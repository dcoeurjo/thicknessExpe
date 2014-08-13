thickness
=========

Travail effectué dans le cadre d'un stage de L3. 

Intitulé du stage : Fonctions d'épaisseur pour l'analyse géométrique de formes.

Contenu du dépôt : Différents tests, pour étudier des fonctions d'épaisseur, ainsi que quelques courbes obtenues et des outils plus généraux ( et plus pratique pour être réutilisés... )

Des Readme ( ou Lisez-moi selon l'humeur ) dans chaque dossier montrent comment se servir des outils ( compilation, prérequis, execution ).

Ce projet utilise quelques outils et librairies, notamment :
- CGAL
- DGtal / DGtalTools
- Binvox
- MAEVA Toolkit

Celles-ci sont rappelées dans les Readme, mais il se peut que des oublis soient survenus.

Pour tout renseignement ou commentaire ( erreur, imprécision... ) : william.aufort@ens-lyon.fr


Descriptif sommaire des tests
==============================

compareNoises : Comparer l'erreur locale causee sur un surface avec différents types de bruits sur les valuers d'épaisseur surfacique (SDF)

compareThickness : Comparer des profils d'épaisseurs surfacique ( SDF ) et volumiques ( avec différentes résolutions )

emd : Calcul de l'Earth Mover Distance entre deux histogrammes.

GlobalError : Calcul de l'erreur globale ( moyenne des erreurs locales )

localErrorSDF : Comparer les profils d'erreurs locales obtenus avec différentes amplitudes de bruits, méthode surfacique

localErrorVol : Idem, mais pour la méthode volumique

Ordonnancement : Contient des outils plus élémentaires qui ont servis dans les différents programmes. 

repartition : Calculer la répartition des rayons à partir des valeurs ( = créer des histogrammes )

results : Des resultats obtenus en utilisant ces outils ( essentiellent des courbes )

