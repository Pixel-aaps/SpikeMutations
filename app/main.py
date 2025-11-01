from fastapi import FastAPI
from pydantic import BaseModel
from typing import Dict
from predictor import predict_mutation

app = FastAPI(title="Spike Protein Mutation Predictor")

# Define expected input dynamically (features must match your model)
class MutationInput(BaseModel):
    features: Dict[str, float]  # Example: {"feature1": 0.23, "feature2": 1.2, ...}

@app.get("/")
def root():
    return {"message": "Spike Protein Mutation Predictor API is running."}

@app.post("/predict")
def predict(data: MutationInput):
    features_dict = data.features
    try:
        prediction = predict_mutation(features_dict)
        return {"prediction": prediction}
    except Exception as e:
        return {"error": str(e)}
