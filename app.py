from flask import Flask, request, jsonify
from flask_cors import CORS
from Bio import Entrez
from googletrans import Translator
import pandas as pd
import os
from datetime import datetime

app = Flask(__name__)
CORS(app)

Entrez.email = "tuemail@gmail.com"  # Usa tu correo válido

@app.route('/buscar', methods=['GET'])
def buscar_articulos():
    categoria = request.args.get('categoria')

    if not categoria:
        return jsonify({"error": "Debes proporcionar una categoría"}), 400

    try:
        # Obtener la fecha actual y establecer el rango desde enero de 2024
        fecha_actual = datetime.today().strftime("%Y/%m/%d")
        fecha_inicio = "2024/01/01"

        # Agregar filtro de fecha a la búsqueda
        query = f'({categoria}) AND ("{fecha_inicio}"[Date - Publication] : "{fecha_actual}"[Date - Publication])'

        handle = Entrez.esearch(db="pubmed", term=query, retmax=20, sort="pub+date")
        record = Entrez.read(handle)
        handle.close()

        article_ids = record.get("IdList", [])
        if not article_ids:
            return jsonify({"error": "No se encontraron artículos desde enero de 2024"}), 404

        resultados = []
        translator = Translator()

        for article_id in article_ids:
            try:
                summary = Entrez.efetch(db="pubmed", id=article_id, rettype="abstract", retmode="xml")
                summary_record = Entrez.read(summary)

                if "PubmedArticle" in summary_record:
                    article = summary_record["PubmedArticle"][0]
                    title = article["MedlineCitation"]["Article"].get("ArticleTitle", "Sin título")
                    abstract = article["MedlineCitation"]["Article"].get("Abstract", {}).get("AbstractText", ["No hay resumen"])[0]
                    pub_date = article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"].get("Year", "Fecha desconocida")
                    link = f"https://pubmed.ncbi.nlm.nih.gov/{article_id}/"

                    # Traducir solo si hay contenido
                    titulo_traducido = translator.translate(title, src="auto", dest="es").text if title else title
                    resumen_traducido = translator.translate(abstract, src="auto", dest="es").text if abstract else abstract

                    resultados.append({
                        "fecha": pub_date,
                        "titulo": titulo_traducido,
                        "resumen": resumen_traducido,
                        "enlace": link
                    })
            except Exception as e:
                print(f"Error procesando artículo {article_id}: {e}")

        return jsonify(resultados)

    except Exception as e:
        return jsonify({"error": f"Error interno en el servidor: {str(e)}"}), 500

if __name__ == '__main__':
    app.run(debug=False, host='0.0.0.0', port=int(os.environ.get("PORT", 5000)))

