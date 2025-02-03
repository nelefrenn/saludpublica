from flask import Flask, request, jsonify, send_file
from Bio import Entrez
from googletrans import Translator
import pandas as pd
import os

app = Flask(__name__)

# Configurar Entrez con un correo electrónico
Entrez.email = "tuemail@ejemplo.com"

@app.route('/buscar', methods=['GET'])
def buscar_articulos():
    categoria = request.args.get('categoria')
    
    # Buscar en PubMed
    handle = Entrez.esearch(db="pubmed", term=categoria, retmax=10, sort="pub+date")
    record = Entrez.read(handle)
    handle.close()
    
    article_ids = record["IdList"]
    handle = Entrez.efetch(db="pubmed", id=",".join(article_ids), rettype="medline", retmode="text")
    data = handle.read()
    handle.close()

    resultados = []
    translator = Translator()
    
    for article_id in article_ids:
        summary = Entrez.efetch(db="pubmed", id=article_id, rettype="abstract", retmode="xml")
        summary_record = Entrez.read(summary)

        if "PubmedArticle" in summary_record:
            article = summary_record["PubmedArticle"][0]
            title = article["MedlineCitation"]["Article"]["ArticleTitle"]
            abstract = article["MedlineCitation"]["Article"].get("Abstract", {}).get("AbstractText", [""])[0]
            pub_date = article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"].get("Year", "No Fecha")
            link = f"https://pubmed.ncbi.nlm.nih.gov/{article_id}/"

            resumen_traducido = translator.translate(abstract, src="auto", dest="es").text

            resultados.append({
                "fecha": pub_date,
                "titulo": title,
                "resumen": resumen_traducido,
                "enlace": link
            })
    
    resultados.sort(key=lambda x: x["fecha"], reverse=True)

    return jsonify(resultados)

@app.route('/descargar_excel')
def descargar_excel():
    data = request.args.get("data")

    if not data:
        return "No hay datos para exportar", 400

    df = pd.DataFrame(data)
    archivo_excel = "resultados_pubmed.xlsx"
    df.to_excel(archivo_excel, index=False)

    return send_file(archivo_excel, as_attachment=True)

if __name__ == '__main__':
    port = int(os.environ.get("PORT", 5000))
    app.run(debug=True, host='0.0.0.0', port=port)
