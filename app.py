from flask import Flask, request, jsonify, send_file
from flask_cors import CORS
from Bio import Entrez
from deep_translator import GoogleTranslator
import pandas as pd
import os
from datetime import datetime

app = Flask(__name__)
CORS(app)

Entrez.email = "tuemail@gmail.com"

@app.route('/descargar_excel', methods=['POST'])
def descargar_excel():
    try:
        datos = request.json.get("datos", [])

        if not datos:
            return jsonify({"error": "No hay datos para exportar"}), 400

        df = pd.DataFrame(datos)
        archivo_excel = "resultados_pubmed.xlsx"
        df.to_excel(archivo_excel, index=False)

        return send_file(archivo_excel, as_attachment=True, mimetype="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")

    except Exception as e:
        return jsonify({"error": f"Error al generar el Excel: {str(e)}"}), 500

if __name__ == '__main__':
    app.run(debug=False, host='0.0.0.0', port=int(os.environ.get("PORT", 5000)))



