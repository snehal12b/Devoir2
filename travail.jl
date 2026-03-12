# ---
# title: Devoir 2
# repository: tpoisot/BIO245-modele
# auteurs:
#    - nom: Bhandari
#      prenom: Snehal
#      matricule: 20279267
#      github: snehal12b
#    - nom: Auteur
#      prenom: Deuxième
#      matricule: XXXXXXXX
#      github: DeuxiAut
# ---

# # Introduction

#La gestion de la végétation sous les lignes électriques à haute tension est considéré comme un défi important.
#Il faut assurer la sécurité des infrastructures et préserver aussi la biodiversité. 
#Une végétation trop dense ou composée d'arbre de grande taille peut interférer avec les lignes électriques, tandis qu'une
#absence complète de végétation peut diminuer la diversité biologique et ainsi favoriser l'érosion des sols.
#Une solution serait de créer des corridors végétalisés principalement composés d'herbes et de buissons de petite taille.
#Cela permet donc de maintenir une couverture végétale et d'éliminer les risques pour les infrastructures.
# # Présentation du modèle

# # Implémentation

# ## Packages nécessaires

import Random
Random.seed!(123456)
using CairoMakie

# ## Une autre section

"""
    foo(x, y)

Cette fonction ne fait rien.
"""
function foo(x, y)
    ## Cette ligne est un commentaire
    return nothing
end

# # Présentation des résultats

# La figure suivante représente des valeurs aléatoires:

hist(randn(100))

# # Discussion



# On peut aussi citer des références dans le document `references.bib`,
# @ermentrout1993cellular -- la bibliographie sera ajoutée automatiquement à la
# fin du document.
