# ---
# title: Devoir 2
# repository: tpoisot/BIO245-modele
# auteurs:
#    - nom: Bhandari
#      prenom: Snehal
#      matricule: 20279267
#      github: snehal12b
#    - nom: Hammel Monzon
#      prenom: Valentina
#      matricule: 20274033
#      github: valentina9000
# ---

# # Introduction

#La gestion de la végétation sous les lignes électriques à haute tension est considéré comme un défi important.
#Il faut assurer la sécurité des infrastructures et préserver aussi la biodiversité. 
#Une végétation trop dense ou composée d'arbre de grande taille peut interférer avec les lignes électriques, tandis qu'une
#absence complète de végétation peut diminuer la diversité biologique et ainsi favoriser l'érosion des sols.
#Une solution serait de créer des corridors végétalisés principalement composés d'herbes et de buissons de petite taille.
#Cela permet donc de maintenir une couverture végétale et d'éliminer les risques pour les infrastructures.
#Ce travail utilise un modèle de transition végétale pour simuler la colonisation et la succession écologique dans un corridor divisé en plusieurs parcelles.
#L'objectif est donc de déterminer une population initiale de plantes et une matrice de transition décrivant les probabilités de changement d'état des parcelles
#qui permettrait d'atteindre un équilibre écologique souhaité. En effet, le mandat indique qu'à l'équilibre, 20% des parcelles doivent être végétalisées, dont 30%
# d'herbes et 70% de buissons. Aussi, pour maintenir une diversité minimale, la variété de buisson la moins abondante doit représenter au moins 30% des parcelles 
#occupées par des buissons. À l'aide de simulations stochastiques et déterministes, nous souhaitons évaluer si les paramètres choisis permettent de respecter ces
#critères dans au moins 80% des simulations.
# # Présentation du modèle

# # Implémentation
# Le modèle à été créer dans le langage Julia afin de simuler l'évolution de la végétation dans un corridor composé de 200 parcelles. Les états que peuvent prendre 
# les parcelles sont les suivant Barren (sol nu), Grass (herbes), Shrub1 (buisson 1) et Shrub2 (buisson 2). 
# Une matrice de transition a été créer poour décrire les probabilités de changement d'états d'une génération à l'autre. Une fonction d'abord vérifie que la matrice 
# de transition est valide en s'assurant que la somme des probabilités de chaque ligne est égale à 1. Sinon, les valeurs sont normalisées. 
# La populaiton initiale est définie par un vecteur contenant le nombre de parcelle dans chaque état. Elle est composée de 160 parcelles nues, 12 parcelles d'herbes, 
# et 14 parcelles pour chacun des deux types de buissons. La quantité d'herbe est volontairement faible pour éviter d'avoir un domiance trop importante d'herbe puisqu'on
# souhaite obtenir 70% de buissons dans la végétation.
# Il y a deux types de simulations, la simulation détéerminsite et la simulation stochastique. La simulation déterministe répresente la tendance moyenne du système 
# alors que la simulation stochastique permet de représenter la variabilité naturelle du système
# Les simulations sont réalisées sur 200 générations. Pour vérifier si les conditions du mandat sont respectées, la simulation stochastique a été répétée 100 fois.
# Après chaque simulation, une fonction vérifie si quatre conditions sont respectées. La première est que le nombre total de parcelle végétalisées doit être environ 
# 40 parcelles avec une tolérance de  + ou - 8 parcelles. La deuxième est que les herbes doivent représenter environ 30 % de la population. La quatrièreme est que 
# les buissons doivent représenter environ 70% de la végétation. La quatrième est que le type de buisson le moins abondant doit représenter au moins 30% des parcelles
# occupées par des buissons. Si les quatres conditions sont respectées alors la simulaiton est réussite. 
# Le pourcentage de simulations qui respecte ces conditions est ensuite calculé pour voir si le modèle atteint l'objectif de 80 % de réussite.
# Finalement, un graphique est créer afin de visualiser l'évolution du nombre de parcelles au cours du temps. Les lignes pales representent les simulation stochastique
# et les lignes épaisses représentent la simulation déterministe. 

# ## Packages nécessaires

import Random
Random.seed!(123456)
using CairoMakie
using Distributions
# ## Code à modifier

# Vérifie que la matrice de transition est valide
function check_transition_matrix!(T)
    for ligne in axes(T, 1)
        if sum(T[ligne, :]) != 1
            @warn "La somme de la ligne $(ligne) n'est pas égale à 1 et a été modifiée"
            T[ligne, :] ./= sum(T[ligne, :])
        end
    end
    return T
end

# Vérifie que les dimensions de la matrice et du vecteur états sont correctes
function check_function_arguments(transitions, states)
    if size(transitions, 1) != size(transitions, 2)
        throw("La matrice de transition n'est pas carrée")
    end

    if size(transitions, 1) != length(states)
        throw("Le nombre d'états ne correspond psa à la matrice de transition")
    end
    return nothing
end

# Simulation stochastique
function _sim_stochastic!(timeseries, transitions, generation)
    for state in axes(timeseries, 1)
        pop_change = rand(Multinomial(timeseries[state, generation], transitions[state, :]))
        timeseries[:, generation+1] .+= pop_change
    end
end

# Simulation déterministe
function _sim_determ!(timeseries, transitions, generation)
    pop_change = (timeseries[:, generation]' * transitions)'
    timeseries[:, generation+1] .= pop_change
end

# Fonction principale de la simulation
function simulation(transitions, states; generations=500, stochastic=false)

    check_transition_matrix!(transitions)
    check_function_arguments(transitions, states)

    _data_type = stochastic ? Int64 : Float32
    timeseries = zeros(_data_type, length(states), generations + 1)
    timeseries[:, 1] = states

    _sim_function! = stochastic ? _sim_stochastic! : _sim_determ!

    for generation in Base.OneTo(generations)
        _sim_function!(timeseries, transitions, generation)
    end

    return timeseries
end

# ## États et population initiale
# States
# Barren, Grass, Shrub1, Shrub2

# Population initiale
s = [150, 0, 25, 25] #200 parcelles et 50 plantées, pas d'herbes initialement parce que l'objectif final est 70% de buissons parmi la végétation, donc si on met les herbes on risque d'en avoir trop à l'équilibre
states = length(s)
patches = sum(s)

# Matrice de transitions
T = zeros(Float64, states, states)
T[1, :] = [150, 12, 6, 6] #Barren, sol nu dominant, mais colonisation possible. Parcelle nue peut devenir herbe, buisson 1, buisson 2
T[2, :] = [25, 95, 10, 10] #Grass, herbes persistent, mais peuvent quand même devenir des buissons
T[3, :] = [10, 8, 110, 12] #Shrub1, buisson 1 persiste
T[4, :] = [10, 8, 12, 110] #Shrub2, buisson 2 persiste, mais il ne domine pas nécessairement le buisson 1, maintient de la diversité
T

states_names = ["Barren", "Grasses", "Shrub1", "Shrub2"]
states_colors = [:grey40, :orange, :teal, :purple]

# Vérification des critères

# Simulations

f = Figure()
ax = Axis(f[1, 1], xlabel="Nb. générations", ylabel="Nb. parcelles")

# Stochastic simulation
for _ in 1:100
    sto_sim = simulation(T, s; stochastic=true, generations=200)
    for i in eachindex(s)
        lines!(ax, sto_sim[i, :], color=states_colors[i], alpha=0.1)
    end
end

# Deterministic simulation
det_sim = simulation(T, s; stochastic=false, generations=200)
for i in eachindex(s)
    lines!(ax, det_sim[i, :], color=states_colors[i], alpha=1, label=states_names[i], linewidth=4)
end

axislegend(ax)
tightlimits!(ax)
current_figure()
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
