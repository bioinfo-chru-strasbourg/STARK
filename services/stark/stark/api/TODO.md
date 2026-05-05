# TODO

Amélioration de l'interface

## Mieux gérer les tasks killed

Lorsque l'action killed est utilisée, il faudrait que la valeur de State le précise. Notamment, le fichier info permet de récupérer le "Exit status".

Exemple de retour de task spooler, le fichier info :
```
Exit status: killed by signal 15
Command: sh -c docker run --rm -v /var/run/docker.sock:/var/run/docker.sock -v /Users/lebechea/STARK/databases:/STARK/databases:ro -v /Users/lebechea/STARK/config/myapps:/STARK/config/myapps:ro -v /Users/lebechea/STARK/config/howard:/STARK/config/howard:ro -v /Users/lebechea/STARK/data:/STARK/data:rw -v /Users/lebechea/STARK/input/runs:/STARK/input/runs:ro -v /Users/lebechea/STARK/input/manifests:/STARK/input/manifests:ro -v /Users/lebechea/STARK/input/pedigree:/STARK/input/pedigree:ro -v /Users/lebechea/STARK/output/repository:/STARK/output/repository:rw -v /Users/lebechea/STARK/output/archives:/STARK/output/archives:rw -v /Users/lebechea/STARK/output/favorites:/STARK/output/favorites:rw -v /Users/lebechea/STARK/services/stark/stark/api:/STARK/services/stark/stark/api -v /STARK/output/results -v /STARK/output/demultiplexing --name STARK.oZGeCmIPejkE.ID-LSubQn5xWW7XvW4U1WpXBYUbZ1zrjpeCNcCw4ZRcS-NAME-External_Analysis_Docker_Command_60  --memory=4G -v /tmp:/tmp --entrypoint=sh  --cpus=2  alpine  -c 'sleep 60 && echo hello from docker' > /STARK/services/stark/stark/api/STARK.oZGeCmIPejkE.ID-LSubQn5xWW7XvW4U1WpXBYUbZ1zrjpeCNcCw4ZRcS-NAME-External_Analysis_Docker_Command_60.output 2>&1 && (echo 'done' > /STARK/services/stark/stark/api/STARK.oZGeCmIPejkE.ID-LSubQn5xWW7XvW4U1WpXBYUbZ1zrjpeCNcCw4ZRcS-NAME-External_Analysis_Docker_Command_60.info && exit 0) || (echo 'failed' > /STARK/services/stark/stark/api/STARK.oZGeCmIPejkE.ID-LSubQn5xWW7XvW4U1WpXBYUbZ1zrjpeCNcCw4ZRcS-NAME-External_Analysis_Docker_Command_60.info && exit 1)
Slots required: 2
Enqueue time: Fri Apr 17 10:58:10 2026
Start time: Fri Apr 17 10:58:10 2026
End time: Fri Apr 17 10:58:23 2026
Time run: 13.195173s
```

## Réorganisation des onlgets/tabs

Il faudrait 4 onglets principaux

### Tasks

L'onglet principal "Tasks" regroupe toutes les tasks sur le cluster. Ce tableau de tasks (tableau déjà implémenté en partie avec le tablean du cluster et le tableau du local) devra :

- lister toutes les tasks de l'ensemble des peers/nodes 
- avoir les colonnes : "Node" "ID" "Queue" "Slots" "State" "Time" "Analysis Name" "Actions"
- être ordonné suivant la State : d'abord "running", puis "queued", puis "finished" (déjà implémenté)
- pouvoir être ordonné suivant les colonnes "Node" "ID" "Queue" "Slots" "State" "Time" "Analysis Name"
- pouvoir filter par les colonnes "Node" "Queue" "State", avec un bouton "reset" (déjà implémenté dans la liste des tasks local)
- pouvoir faire une recherche sur la colonne "Analysis Name"

### Cluster

Un résumé de l'état du cluster en terme de ressources (déjà implémenté).

### Launch analysis

Permet de lancer une analyse (déjà implémenté dans l'onglet "Local").

- Un onglet "Run" permettra de lancer un Run par son nom (déjà implémenté)
- Un onglet "Advanced" permettra de lancer une analyse par un JSON (déjà implémenté)

Une fois l'analyse lancée, un récapitulatif clair devra indiquer ce qui a été lancé :

- l'identifiant de l'analyse (exemple "STARK.arMryjKRfZeb.ID-5KQNEgo5nb1SIngtx9zO97HpRykEQt5UONoDYFgUI-NAME-External_Analysis_Docker_Command_60E2")
- sur quel peer/node
- sur quel queue
- quel ID
- combien de slots (sur combien configuré)

### Archives

Cet onglet permet de récupérer toutes les tasks qui auront été lancées, en se basant uniquement sur les fichiers générés (info, json, output).
Un tableau sera généré, de la même manière que pour le tableau des tasks
