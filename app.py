from flask import Flask, request, jsonify, send_file
from flask_cors import CORS
from Bio import Entrez
from googletrans import Translator
import pandas as pd
import os
from datetime import datetime

app = Flask(__name__)
CORS(app)

Entrez.email = "nelefren@gmail.com"

@app.route('/buscar', methods=['GET'])
def buscar_articulos():
    categoria = request.args.get('categoria')
    fecha_actual = datetime.today().strftime("%Y/%m/%d")
    fecha_inicio = "2024/01/01"
    query = f'({categoria}) AND ("{fecha_inicio}"[Date - Publication] : "{fecha_actual}"[Date - Publication])'

    handle = Entrez.esearch(db="pubmed", term=query, retmax=20, sort="pub+date")
    record = Entrez.read(handle)
    handle.close()

    article_ids = record.get("IdList", [])
    if not article_ids:
        return jsonify({"error": "No se encontraron artículos"}), 404

    return jsonify(article_ids)

@app.route('/descargar_excel', methods=['POST'])
def descargar_excel():
    datos = request.json.get("datos", [])
    df = pd.DataFrame(datos)
    archivo_excel = "resultados_pubmed.xlsx"
    df.to_excel(archivo_excel, index=False)
    return send_file(archivo_excel, as_attachment=True)

if __name__ == '__main__':
    app.run(debug=False, host='0.0.0.0', port=int(os.environ.get("PORT", 5000)))
