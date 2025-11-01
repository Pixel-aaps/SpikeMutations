# train_model.py
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.preprocessing import StandardScaler, LabelEncoder, OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix
import joblib
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset

df = pd.read_csv("final_dataset_ready.csv")

TARGET = "label"  
X = df.drop(columns=[TARGET])
y = df[TARGET]

if y.dtype == 'object':
    le = LabelEncoder()
    y = le.fit_transform(y)
    joblib.dump(le, "label_encoder.pkl")

numeric_features = X.select_dtypes(include=[np.number]).columns.tolist()
categorical_features = X.select_dtypes(exclude=[np.number]).columns.tolist()

print("Numeric features:", numeric_features)
print("Categorical features:", categorical_features)

numeric_transformer = Pipeline(steps=[
    ("imputer", SimpleImputer(strategy="mean")),
    ("scaler", StandardScaler())
])

categorical_transformer = Pipeline(steps=[
    ("imputer", SimpleImputer(strategy="most_frequent")),
    ("onehot", OneHotEncoder(handle_unknown="ignore"))
])

preprocessor = ColumnTransformer(
    transformers=[
        ("num", numeric_transformer, numeric_features),
        ("cat", categorical_transformer, categorical_features),
    ]
)

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, stratify=y, random_state=42
)

rf_model = RandomForestClassifier(class_weight="balanced", random_state=42)

pipe = Pipeline(steps=[("preprocessor", preprocessor),
                       ("classifier", rf_model)])

param_grid = {
    "classifier__n_estimators": [100, 200],
    "classifier__max_depth": [None, 10, 20],
}

grid = GridSearchCV(pipe, param_grid, cv=3, n_jobs=-1, scoring="f1_macro")
grid.fit(X_train, y_train)

best_model = grid.best_estimator_

y_pred = best_model.predict(X_test)
print("Classification Report:\n", classification_report(y_test, y_pred))
print("Confusion Matrix:\n", confusion_matrix(y_test, y_pred))

joblib.dump(best_model, "mutation_model.pkl")
print("RandomForest model saved as mutation_model.pkl")

# Transform dataset to numeric only for DL
X_proc = preprocessor.fit_transform(X)
X_train_dl, X_test_dl, y_train_dl, y_test_dl = train_test_split(
    X_proc, y, test_size=0.2, stratify=y, random_state=42
)

# Convert to tensors
X_train_tensor = torch.tensor(
    X_train_dl.toarray() if hasattr(X_train_dl, "toarray") else X_train_dl,
    dtype=torch.float32
)
X_test_tensor = torch.tensor(
    X_test_dl.toarray() if hasattr(X_test_dl, "toarray") else X_test_dl,
    dtype=torch.float32
)
y_train_tensor = torch.tensor(y_train_dl.values, dtype=torch.long)
y_test_tensor  = torch.tensor(y_test_dl.values, dtype=torch.long)


train_dataset = TensorDataset(X_train_tensor, y_train_tensor)
test_dataset = TensorDataset(X_test_tensor, y_test_tensor)
train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)
test_loader = DataLoader(test_dataset, batch_size=32, shuffle=False)

# Define simple MLP
class MLPClassifier(nn.Module):
    def __init__(self, input_dim, hidden_dim=128, output_dim=len(np.unique(y))):
        super(MLPClassifier, self).__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(hidden_dim, output_dim)
        )
    def forward(self, x):
        return self.net(x)

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = MLPClassifier(X_train_tensor.shape[1]).to(device)
criterion = nn.CrossEntropyLoss()
optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)

# Training loop
EPOCHS = 20
for epoch in range(EPOCHS):
    model.train()
    for X_batch, y_batch in train_loader:
        X_batch, y_batch = X_batch.to(device), y_batch.to(device)
        optimizer.zero_grad()
        outputs = model(X_batch)
        loss = criterion(outputs, y_batch)
        loss.backward()
        optimizer.step()
    print(f"Epoch {epoch+1}/{EPOCHS}, Loss: {loss.item():.4f}")

model.eval()
y_pred_dl = []
with torch.no_grad():
    for X_batch, _ in test_loader:
        X_batch = X_batch.to(device)
        outputs = model(X_batch)
        y_pred_dl.extend(torch.argmax(outputs, dim=1).cpu().numpy())

print("PyTorch model evaluation")
print(classification_report(y_test_dl, y_pred_dl))