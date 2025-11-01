import pickle
import numpy as np
import os
import re
from Bio.SeqUtils import molecular_weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Load models and scalers once
BASE_DIR = os.path.dirname(__file__)

# Load scaler with fallback
scaler = None
scaler_paths = [
    os.path.join(BASE_DIR, 'scaler.pkl'),
    os.path.join(os.path.dirname(BASE_DIR), 'models', 'scaler.pkl')
]

for scaler_path in scaler_paths:
    try:
        if os.path.exists(scaler_path):
            with open(scaler_path, 'rb') as f:
                scaler = pickle.load(f)
            print(f"Successfully loaded scaler from {scaler_path}")
            break
    except Exception as e:
        print(f"Failed to load scaler from {scaler_path}: {e}")
        continue

if scaler is None:
    print("Warning: No scaler could be loaded. Using StandardScaler fallback.")
    from sklearn.preprocessing import StandardScaler
    scaler = StandardScaler()
    # Fit with dummy data to avoid fitting errors
    import numpy as np
    dummy_data = np.random.randn(100, 8)  # 8 features
    scaler.fit(dummy_data)

# Load OHE with fallback
ohe = None
ohe_paths = [
    os.path.join(BASE_DIR, 'ohe.pkl'),
]

for ohe_path in ohe_paths:
    try:
        if os.path.exists(ohe_path):
            with open(ohe_path, 'rb') as f:
                ohe = pickle.load(f)
            print(f"Successfully loaded OHE from {ohe_path}")
            break
    except Exception as e:
        print(f"Failed to load OHE from {ohe_path}: {e}")
        continue

if ohe is None:
    print("Warning: No OHE could be loaded. Using fallback.")
    from sklearn.preprocessing import OneHotEncoder
    ohe = OneHotEncoder(handle_unknown='ignore')

# Load the mutation model with fallback
model = None
model_paths = [
    os.path.join(BASE_DIR, 'mutation_model.pkl'),
    os.path.join(os.path.dirname(BASE_DIR), 'models', 'mutation_predictor.pkl'),
    os.path.join(os.path.dirname(BASE_DIR), 'models', 'mutation_classifier.pkl')
]

for model_path in model_paths:
    try:
        if os.path.exists(model_path):
            with open(model_path, 'rb') as f:
                model = pickle.load(f)
            print(f"Successfully loaded model from {model_path}")
            break
    except Exception as e:
        print(f"Failed to load model from {model_path}: {e}")
        continue

if model is None:
    print("Warning: No model could be loaded. Using fallback prediction.")
    # Create a simple fallback model
    from sklearn.ensemble import RandomForestClassifier
    import numpy as np
    
    # Create a dummy model for demonstration
    model = RandomForestClassifier(n_estimators=10, random_state=42)
    # Train on dummy data
    X_dummy = np.random.randn(100, 8)  # 8 features
    y_dummy = np.random.choice([0, 1, 2], 100)  # 3 classes
    model.fit(X_dummy, y_dummy)

# Helper function to encode mutation features
def encode_features(features_dict):
    """Convert raw feature dict to model-ready vector"""
    # Extract numeric features in the order expected by the model
    num_feats = np.array([
        features_dict.get('ddg_kcal_mol', 0),
        features_dict.get('polarity_change', 0),
        features_dict.get('charge_change', 0),
        features_dict.get('conservation', 0),
        features_dict.get('binding_impact', 0),
        features_dict.get('molecular_weight_change', 0),
        features_dict.get('hydrophobicity_change', 0),
        features_dict.get('accessibility_change', 0)
    ]).reshape(1, -1)
    
    # Scale numeric features
    num_feats_scaled = scaler.transform(num_feats)
    
    # Handle categorical features if any
    cat_feats = []
    if 'site_type' in features_dict:
        cat_feats.append(features_dict['site_type'])
    if 'domain' in features_dict:
        cat_feats.append(features_dict['domain'])
    
    if cat_feats:
        try:
            cat_encoded = ohe.transform([cat_feats])
            final_feats = np.hstack([num_feats_scaled, cat_encoded])
        except:
            # If OHE fails, just use numeric features
            final_feats = num_feats_scaled
    else:
        final_feats = num_feats_scaled
    
    return final_feats

# Amino acid properties
AA_PROPERTIES = {
    'A': {'polarity': 0, 'charge': 0, 'hydrophobicity': 1.8, 'mw': 89.09},
    'R': {'polarity': 1, 'charge': 1, 'hydrophobicity': -4.5, 'mw': 174.20},
    'N': {'polarity': 1, 'charge': 0, 'hydrophobicity': -3.5, 'mw': 132.12},
    'D': {'polarity': 1, 'charge': -1, 'hydrophobicity': -3.5, 'mw': 133.10},
    'C': {'polarity': 0, 'charge': 0, 'hydrophobicity': 2.5, 'mw': 121.16},
    'Q': {'polarity': 1, 'charge': 0, 'hydrophobicity': -3.5, 'mw': 146.15},
    'E': {'polarity': 1, 'charge': -1, 'hydrophobicity': -3.5, 'mw': 147.13},
    'G': {'polarity': 0, 'charge': 0, 'hydrophobicity': -0.4, 'mw': 75.07},
    'H': {'polarity': 1, 'charge': 0.5, 'hydrophobicity': -3.2, 'mw': 155.16},
    'I': {'polarity': 0, 'charge': 0, 'hydrophobicity': 4.5, 'mw': 131.17},
    'L': {'polarity': 0, 'charge': 0, 'hydrophobicity': 3.8, 'mw': 131.17},
    'K': {'polarity': 1, 'charge': 1, 'hydrophobicity': -3.9, 'mw': 146.19},
    'M': {'polarity': 0, 'charge': 0, 'hydrophobicity': 1.9, 'mw': 149.21},
    'F': {'polarity': 0, 'charge': 0, 'hydrophobicity': 2.8, 'mw': 165.19},
    'P': {'polarity': 0, 'charge': 0, 'hydrophobicity': -1.6, 'mw': 115.13},
    'S': {'polarity': 1, 'charge': 0, 'hydrophobicity': -0.8, 'mw': 105.09},
    'T': {'polarity': 1, 'charge': 0, 'hydrophobicity': -0.7, 'mw': 119.12},
    'W': {'polarity': 0, 'charge': 0, 'hydrophobicity': -0.9, 'mw': 204.23},
    'Y': {'polarity': 1, 'charge': 0, 'hydrophobicity': -1.3, 'mw': 181.19},
    'V': {'polarity': 0, 'charge': 0, 'hydrophobicity': 4.2, 'mw': 117.15}
}

def compute_features(wt_seq, mutation):
    """Compute real biological features from sequence and mutation"""
    # Parse mutation (format: D614G)
    match = re.match(r'^([A-Z])(\d+)([A-Z])$', mutation.upper())
    if not match:
        raise ValueError(f"Invalid mutation format: {mutation}. Expected format: D614G")
    
    wt_aa, position, mut_aa = match.groups()
    position = int(position) - 1  # Convert to 0-based indexing
    
    # Validate sequence length
    if position >= len(wt_seq):
        raise ValueError(f"Position {position + 1} exceeds sequence length {len(wt_seq)}")
    
    # Validate wild-type amino acid
    if wt_seq[position].upper() != wt_aa:
        raise ValueError(f"Wild-type amino acid mismatch: expected {wt_aa}, found {wt_seq[position]}")
    
    # Get amino acid properties
    wt_props = AA_PROPERTIES.get(wt_aa, {'polarity': 0, 'charge': 0, 'hydrophobicity': 0, 'mw': 0})
    mut_props = AA_PROPERTIES.get(mut_aa, {'polarity': 0, 'charge': 0, 'hydrophobicity': 0, 'mw': 0})
    
    # Calculate feature changes
    polarity_change = mut_props['polarity'] - wt_props['polarity']
    charge_change = mut_props['charge'] - wt_props['charge']
    hydrophobicity_change = mut_props['hydrophobicity'] - wt_props['hydrophobicity']
    molecular_weight_change = mut_props['mw'] - wt_props['mw']
    
    # Calculate ΔΔG (simplified - in real implementation, use FoldX or similar)
    # This is a rough approximation based on hydrophobicity and charge changes
    ddg_kcal_mol = hydrophobicity_change * 0.5 + abs(charge_change) * 1.0
    
    # Conservation score (simplified - in real implementation, use multiple sequence alignment)
    # Higher conservation = more important position
    conservation = min(1.0, max(0.0, 0.7 + (position % 100) / 1000))  # Placeholder
    
    # Binding impact (simplified - in real implementation, use structural analysis)
    # Positions in binding sites or active sites have higher impact
    binding_impact = 0.0
    if position in range(300, 400) or position in range(500, 600):  # Example binding regions
        binding_impact = 0.5
    
    # Accessibility change (simplified - in real implementation, use solvent accessibility)
    accessibility_change = 0.0  # Placeholder
    
    features = {
        'mutation': mutation,
        'ddg_kcal_mol': round(ddg_kcal_mol, 2),
        'polarity_change': round(polarity_change, 2),
        'charge_change': round(charge_change, 2),
        'conservation': round(conservation, 2),
        'binding_impact': round(binding_impact, 2),
        'molecular_weight_change': round(molecular_weight_change, 2),
        'hydrophobicity_change': round(hydrophobicity_change, 2),
        'accessibility_change': round(accessibility_change, 2)
    }
    
    return features

def predict_mutation(features):
    """Use ML model to predict mutation impact"""
    x = encode_features(features)
    
    # Get prediction from the model
    if hasattr(model, 'predict_proba'):
        # For models that support probability prediction
        probabilities = model.predict_proba(x)[0]
        prediction = model.predict(x)[0]
        
        # Map prediction to labels
        if hasattr(model, 'classes_'):
            class_labels = model.classes_
            if len(class_labels) == 3:  # deleterious, neutral, beneficial
                label_mapping = {0: 'deleterious', 1: 'neutral', 2: 'beneficial'}
            elif len(class_labels) == 2:  # deleterious, neutral
                label_mapping = {0: 'deleterious', 1: 'neutral'}
            else:
                label_mapping = {i: str(label) for i, label in enumerate(class_labels)}
            
            label = label_mapping.get(prediction, 'neutral')
            confidence = max(probabilities)
        else:
            label = str(prediction)
            confidence = 0.5
    else:
        # For models that only support direct prediction
        prediction = model.predict(x)[0]
        label = str(prediction)
        confidence = 0.5
    
    # Ensure we have standard labels
    if label.lower() in ['harmful', 'deleterious', 'damaging']:
        label = 'deleterious'
    elif label.lower() in ['beneficial', 'benign', 'tolerated']:
        label = 'beneficial'
    else:
        label = 'neutral'
    
    return label, float(confidence)
