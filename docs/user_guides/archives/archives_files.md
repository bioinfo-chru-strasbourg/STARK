# Archives files


```python exec="on"
import plotly.offline as offline
import plotly.graph_objs as go
import yaml 

# data = go.Scatter(x=[1, 2, 3, 4], y=[4, 3, 2, 1])
# layout = go.Layout(title="hello world")
# figure = go.Figure(data=[data], layout=layout)
# print(offline.plot(figure, output_type='div'))


# List all files in pattern /STARK/output/repository/*/*/*/STARComplete.txt

import glob
import os
import datetime

print("TEST")
print("")

# Le chemin du répertoire courant, relatif au dossier racine du projet (où mkdocs.yml se trouve)
# Ce chemin doit être ajusté pour chaque fichier.
current_dir = '/tool/docs/user_guides/archives'
current_relative_dir = 'user_guides/archives'

try:
    # On liste les fichiers dans le répertoire spécifié
    for f in sorted(os.listdir(current_dir)):
        # On évite d'afficher le fichier markdown lui-même
        #if f != 'archives_files.md' and not os.path.isdir(os.path.join(current_dir,f)):
        if not f.endswith(".md") and not os.path.isdir(os.path.join(current_dir,f)):
            #     f_path = f.replace(".md", "")
            # else:
            #     f_path = f
            print(f"- [{f}](/{current_relative_dir}/{f})")
except FileNotFoundError:
    print(f"Erreur : Le répertoire '{current_dir}' n'a pas été trouvé.")
    print("Veuillez vérifier le chemin dans le script.")

```

