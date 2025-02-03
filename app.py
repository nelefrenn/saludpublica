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
    
    if not categoria:
        return jsonify({"error": "Debe proporcionar una categoría de búsqueda"}), 400

    translator = Translator()
    
    try:
        # Traducir el término de búsqueda al inglés
        categoria_en = translator.translate(categoria, src="es", dest="en").text
    except Exception as e:
        return jsonify({"error": f"Error en la traducción: {str(e)}"}), 500
    
    # Buscar en PubMed usando el término traducido
    handle = Entrez.esearch(db="pubmed", term=categoria_en, retmax=10, sort="pub+date")
    record = Entrez.read(handle)
    handle.close()
    
    article_ids = record.get("IdList", [])
    
    if not article_ids:
        return jsonify({"message": "No se encontraron artículos en PubMed"}), 404

    resultados = []
    
    for article_id in article_ids:
        summary = Entrez.efetch(db="pubmed", id=article_id, rettype="abstract", retmode="xml")
        summary_record = Entrez.read(summary)
        summary.close()

        if "PubmedArticle" in summary_record:
            article = summary_record["PubmedArticle"][0]
            title = article["MedlineCitation"]["Article"]["ArticleTitle"]
            abstract = article["MedlineCitation"]["Article"].get("Abstract", {}).get("AbstractText", [""])[0]
            pub_date = article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"].get("Year", "No Fecha")
            link = f"https://pubmed.ncbi.nlm.nih.gov/{article_id}/"

            # Traducir el resumen al español
            try:
                resumen_traducido = translator.translate(abstract, src="en", dest="es").text
            except Exception:
                resumen_traducido = "Traducción no disponible"

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
