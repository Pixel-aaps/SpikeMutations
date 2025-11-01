from flask import Flask, render_template, request, jsonify, send_from_directory
from predictor import compute_features, predict_mutation
from drug_predictor import predict_drug_binding, get_drug_binding_summary
import os
import logging
from logging.handlers import RotatingFileHandler

app = Flask(__name__)
# Resolve absolute paths so files serve correctly regardless of CWD
BASE_DIR = os.path.dirname(__file__)
ROOT_DIR = os.path.dirname(BASE_DIR)
STRUCT_DIR = os.path.join(BASE_DIR, 'structures')
WT_STRUCT_DIR = os.path.join(ROOT_DIR, 'data', 'external', 'Spike3Dstructure')
MUT_STRUCT_DIR = os.path.join(ROOT_DIR, 'data', 'external', 'Mutant3DStructures')

# -------------------------
# Logging configuration
# -------------------------
LOG_DIR = os.path.join(BASE_DIR, 'logs')
os.makedirs(LOG_DIR, exist_ok=True)

file_handler = RotatingFileHandler(os.path.join(LOG_DIR, 'app.log'), maxBytes=1_000_000, backupCount=3)
file_handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(name)s %(message)s'))
file_handler.setLevel(logging.INFO)

if not any(isinstance(h, RotatingFileHandler) for h in app.logger.handlers):
    app.logger.addHandler(file_handler)

app.logger.setLevel(logging.INFO)
app.logger.info('Application starting...')
app.logger.info(f"STRUCT_DIR: {STRUCT_DIR}")
app.logger.info(f"WT_STRUCT_DIR: {WT_STRUCT_DIR}")
app.logger.info(f"MUT_STRUCT_DIR: {MUT_STRUCT_DIR}")

def generate_explanation(features, ml_label, ml_score):
    """Generate a detailed explanation of the mutation prediction"""
    mutation = features['mutation']
    ddg = features['ddg_kcal_mol']
    polarity_change = features['polarity_change']
    charge_change = features['charge_change']
    conservation = features['conservation']
    
    explanation = f"Mutation {mutation} Analysis:\n\n"
    
    # Stability analysis
    if ddg > 1.0:
        explanation += f"• Stability Impact: Destabilizing (ΔΔG = {ddg:.2f} kcal/mol)\n"
    elif ddg < -1.0:
        explanation += f"• Stability Impact: Stabilizing (ΔΔG = {ddg:.2f} kcal/mol)\n"
    else:
        explanation += f"• Stability Impact: Neutral (ΔΔG = {ddg:.2f} kcal/mol)\n"
    
    # Chemical properties
    if polarity_change != 0:
        explanation += f"• Polarity Change: {polarity_change:+.1f}\n"
    if charge_change != 0:
        explanation += f"• Charge Change: {charge_change:+.1f}\n"
    
    # Conservation
    if conservation > 0.8:
        explanation += f"• Conservation: High ({conservation:.2f}) - functionally important\n"
    elif conservation > 0.5:
        explanation += f"• Conservation: Moderate ({conservation:.2f})\n"
    else:
        explanation += f"• Conservation: Low ({conservation:.2f})\n"
    
    # ML prediction
    explanation += f"\nML Prediction: {ml_label.title()} (confidence: {ml_score:.2f})\n"
    
    if ml_label == 'deleterious':
        explanation += "This mutation is predicted to be harmful to protein function."
    elif ml_label == 'beneficial':
        explanation += "This mutation is predicted to be beneficial to protein function."
    else:
        explanation += "This mutation is predicted to have neutral effects on protein function."
    
    return explanation

@app.route('/structures/<filename>')
def serve_structure(filename):
    """Serve structure files"""
    path = os.path.join(STRUCT_DIR, filename)
    if not os.path.exists(path):
        app.logger.warning(f"Requested structure not found: {path}")
    else:
        app.logger.info(f"Serving structure: {path}")
    return send_from_directory(STRUCT_DIR, filename)

@app.route('/structures/wt/<filename>')
def serve_wt_structure(filename):
    """Serve wild-type structure files from data/external/Spike3Dstructure"""
    path = os.path.join(WT_STRUCT_DIR, filename)
    if not os.path.exists(path):
        app.logger.warning(f"Requested WT structure not found: {path}")
    else:
        app.logger.info(f"Serving WT structure: {path}")
    return send_from_directory(WT_STRUCT_DIR, filename)

@app.route('/structures/mutant/<path:filename>')
def serve_mutant_structure(filename):
    """Serve curated mutant structure files (supports nested subfolders)."""
    path = os.path.join(MUT_STRUCT_DIR, filename)
    if not os.path.exists(path):
        app.logger.warning(f"Requested mutant structure not found: {path}")
    else:
        app.logger.info(f"Serving mutant structure: {path}")
    directory = os.path.dirname(path)
    basename = os.path.basename(path)
    return send_from_directory(directory, basename)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/api/analyze', methods=['POST'])
def analyze():
    data = request.get_json()
    wt_seq = data.get('wt_sequence')
    mutation = data.get('mutation')

    if not wt_seq or not mutation:
        return jsonify({'error': 'Missing WT sequence or mutation'}), 400

    try:
        # 1. Compute features
        features = compute_features(wt_seq, mutation)

        # 2. Predict mutation impact
        ml_label, ml_score = predict_mutation(features)
        features['ml_prediction'] = ml_label
        features['ml_score'] = ml_score

        # 3. Provide AI explanation
        explanation = generate_explanation(features, ml_label, ml_score)

        # 4. Predict drug binding
        drug_predictions = get_drug_binding_summary(mutation)
        # If error, just return as is
        if "error" in drug_predictions:
            return jsonify({'error': drug_predictions["error"]}), 500

        # Convert all_predictions dict to a list of dicts with status
        all_preds = []
        for drug, score in drug_predictions.get("all_predictions", {}).items():
            if score > 0.7:
                status = 'High'
            elif 0.5 <= score <= 0.7:
                status = 'Mid'
            else:
                status = 'Low'
            all_preds.append({
                "drug": drug,
                "binding_score": score,
                "status": status
            })
        drug_predictions["all_predictions"] = all_preds

        # 5. Serve PDB structures
        # Default WT structure: use the first AlphaFold model from external data if present
        default_wt = "fold_2025_07_21_20_47_model_0.cif"
        wt_external_path = os.path.join(WT_STRUCT_DIR, default_wt)
        if os.path.exists(wt_external_path):
            wildtype_pdb = f"/structures/wt/{default_wt}"
            app.logger.info(f"Using WT external structure: {wt_external_path}")
        else:
            wildtype_pdb = f"/structures/spike_model_0.cif"  # fallback if not present
            app.logger.info("External WT structure not found, using fallback spike_model_0.cif from /structures")
        # Try curated mutant first
        mutant_pdb = None
        try:
            folder = os.path.join(MUT_STRUCT_DIR, f"{mutation}_Mutant")
            if os.path.isdir(folder):
                picks = [f for f in os.listdir(folder) if f.endswith('_model_0.cif')]
                if not picks:
                    picks = [f for f in os.listdir(folder) if f.lower().endswith('.cif')]
                if picks:
                    mutant_pdb = f"/structures/mutant/{os.path.basename(folder)}/{picks[0]}"
                    app.logger.info(f"Curated mutant found for {mutation}: {mutant_pdb}")
        except Exception as _e:
            app.logger.warning(f"Curated mutant lookup failed for {mutation}: {_e}")

        if mutant_pdb is None:
            # bundled local path
            bundled = os.path.join(STRUCT_DIR, f"spike_{mutation}.cif")
            if os.path.exists(bundled):
                mutant_pdb = f"/structures/spike_{mutation}.cif"
            else:
                mutant_pdb = wildtype_pdb

        return jsonify({
            'features': features,
            'explanation': explanation,
            'drug_predictions': drug_predictions,
            'wildtype_pdb': wildtype_pdb,
            'mutant_pdb': mutant_pdb
        })
    except Exception as e:
        app.logger.exception(f"/api/analyze failed: {e}")
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    app.logger.info('Starting Flask dev server on http://127.0.0.1:5000')
    app.run(debug=True)
