from flask import Flask, request, jsonify, send_file
from flask_cors import CORS
from Bio import Entrez
from googletrans import Translator
import pandas as pd
import os

app = Flask(__name__)
CORS(app, resources={r"/*": {"origins": "*"}})

Entrez.email = "nelefren@gmail.com"  # Pon tu correo real aquí

@app.route('/buscar', methods=['GET'])
def buscar_articulos():
    categoria = request.args.get('categoria')

    if not categoria:
        return jsonify({"error": "Debes proporcionar una categoría"}), 400

    try:
        handle = Entrez.esearch(db="pubmed", term=categoria, retmax=10, sort="pub+date")
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

                # Obtener título y resumen
                title = article["MedlineCitation"]["Article"].get("ArticleTitle", "Sin título")
                abstract = article["MedlineCitation"]["Article"].get("Abstract", {}).get("AbstractText", ["No hay resumen"])[0]
                pub_date = article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"].get("Year", "Fecha desconocida")
                link = f"https://pubmed.ncbi.nlm.nih.gov/{article_id}/"

                # Traducir título y resumen
                try:
                    titulo_traducido = translator.translate(title, src="auto", dest="es").text
                except Exception:
                    titulo_traducido = title  # Si falla, usa el original

                try:
                    resumen_traducido = translator.translate(abstract, src="auto", dest="es").text
                except Exception:
                    resumen_traducido = abstract  # Si falla, usa el original

                resultados.append({
                    "fecha": pub_date,
                    "titulo": titulo_traducido,
                    "resumen": resumen_traducido,
                    "enlace": link
                })

        resultados.sort(key=lambda x: x["fecha"], reverse=True)
        return jsonify(resultados)

    except Exception as e:
        return jsonify({"error": f"Error interno: {str(e)}"}), 500

if __name__ == '__main__':
    port = int(os.environ.get("PORT", 5000))
    app.run(debug=True, host='0.0.0.0', port=port)


