from flask import Flask, request, jsonify, send_file
from flask_cors import CORS  # Importar CORS
from Bio import Entrez
from googletrans import Translator
import pandas as pd
import os

app = Flask(__name__)

# Habilitar CORS correctamente
CORS(app)  # Esto permitirá que cualquier dominio acceda a la API

Entrez.email = "tuemail@gmail.com"  # Usa un correo válido

@app.route('/buscar', methods=['GET'])
def buscar_articulos():
    categoria = request.args.get('categoria')

    if not categoria:
        return jsonify({"error": "Debes proporcionar una categoría"}), 400

    try:
        handle = Entrez.esearch(db="pubmed", term=categoria, retmax=50, sort="pub+date")
        record = Entrez.read(handle)
        handle.close()

        article_ids = record.get("IdList", [])
        if not article_ids:
            return jsonify({"error": "No se encontraron artículos"}), 404

        resultados = []
        translator = Translator()

        for article_id in article_ids:
            summary = Entrez.efetch(db="pubmed", id=article_id, rettype="abstract", retmode="xml")
            summary_record = Entrez.read(summary)

            if "PubmedArticle" in summary_record:
                article = summary_record["PubmedArticle"][0]
                title = article["MedlineCitation"]["Article"].get("ArticleTitle", "Sin título")
                abstract = article["MedlineCitation"]["Article"].get("Abstract", {}).get("AbstractText", ["No hay resumen"])[0]
                pub_date = article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"].get("Year", "Fecha desconocida")
                link = f"https://pubmed.ncbi.nlm.nih.gov/{article_id}/"

                # Agregar manejo de errores en la traducción
                try:
                    titulo_traducido = translator.translate(title, src="auto", dest="es").text
                except Exception:
                    titulo_traducido = title  # En caso de error, usa el original

                try:
                    resumen_traducido = translator.translate(abstract, src="auto", dest="es").text
                except Exception:
                    resumen_traducido = abstract  # En caso de error, usa el original

                resultados.append({
                    "fecha": pub_date,
                    "titulo": titulo_traducido,
                    "resumen": resumen_traducido,
                    "enlace": link
                })

        return jsonify(resultados)

    except Exception as e:
        return jsonify({"error": f"Error interno en el servidor: {str(e)}"}), 500

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=int(os.environ.get("PORT", 5000)))

