from flask import Flask, request, jsonify, send_file
from flask_cors import CORS
from Bio import Entrez
from deep_translator import GoogleTranslator
import pandas as pd
import os
from datetime import datetime

app = Flask(__name__)
CORS(app, resources={r"/*": {"origins": "*"}})  # Habilita CORS para todas las rutas

Entrez.email = "nelefren@gmail.com"  # Asegúrate de usar un correo válido

@app.route('/')
def home():
    return jsonify({"mensaje": "API funcionando correctamente en Render"}), 200

@app.route('/buscar', methods=['GET'])
def buscar_articulos():
    categoria = request.args.get('categoria')
    if not categoria:
        return jsonify({"error": "Debes proporcionar una categoría"}), 400

    try:
        # Filtrar artículos desde enero de 2024 hasta la fecha actual
        fecha_actual = datetime.today().strftime("%Y/%m/%d")
        fecha_inicio = "2024/01/01"
        query = f'({categoria}) AND ("{fecha_inicio}"[Date - Publication] : "{fecha_actual}"[Date - Publication])'

        # Obtener IDs de artículos en PubMed
        handle = Entrez.esearch(db="pubmed", term=query, retmax=10, sort="pub+date")
        record = Entrez.read(handle)
        handle.close()

        article_ids = record.get("IdList", [])
        if not article_ids:
            return jsonify({"error": "No se encontraron artículos desde enero de 2024"}), 404

        resultados = []
        translator = GoogleTranslator(source="auto", target="es")

        # Obtener los detalles de los artículos usando Entrez.efetch()
        handle = Entrez.efetch(db="pubmed", id=",".join(article_ids), rettype="xml", retmode="text")
        article_data = Entrez.read(handle)
        handle.close()

        # Extraer información de cada artículo
        for article in article_data["PubmedArticle"]:
            try:
                title = article["MedlineCitation"]["Article"].get("ArticleTitle", "Sin título")
                abstract_list = article["MedlineCitation"]["Article"].get("Abstract", {}).get("AbstractText", [])
                abstract = abstract_list[0] if abstract_list else "No hay resumen"
                pub_date = article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"].get("Year", "Fecha desconocida")
                article_id = article["MedlineCitation"]["PMID"]
                link = f"https://pubmed.ncbi.nlm.nih.gov/{article_id}/"

                # Traducir título y resumen
                titulo_traducido = translator.translate(title) if title else title
                resumen_traducido = translator.translate(abstract) if abstract else abstract

                resultados.append({
                    "fecha": pub_date,
                    "titulo": titulo_traducido,
                    "resumen": resumen_traducido,
                    "enlace": link
                })
            except Exception as e:
                print(f"Error procesando artículo: {e}")

        return jsonify(resultados)

    except Exception as e:
        return jsonify({"error": f"Error interno en el servidor: {str(e)}"}), 500

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




