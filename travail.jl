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

# Le corridor étudié est représenté comme un ensemble de 200 parcelles qui peuvent se trouver dans différents états de végétation. Quatre états sont considérés : Barren (sol nu), Grass (herbes), Shrub1 (buisson 1) et Shrub2 (buisson 2).
# Les parcelles évoluent d'un état à l'autre au fil du temps selon une matrice de transition qui décrit les probabilités de colonisation, persistance ou de remplacement entre les différents types de végétation.
# Le modèle est un modèle de transition markovien où l'état d'une parcelle à la génération suivante dépend juste de son état actuel et des probabilités associées dans la matrice de transition. 
# Les transitions permettent la colonisation des parcelles nues par des herbes (Grass) ou des buissons (Shrub1 et 2), la persistance des végétaux, et les conversions entre les différents types de végétation.
# Ainsi, pour simuler l'aménagement du corridor sous une ligne électique, une population initiale comportement principalement des parcelles nues est utilisée, avec une plantation maximale de 50 parcelles sous forme de buissons.
# Les herbes ne sont pas plantées initialement, car elles devraient coloniser naturellement les parcelles libres. Donc, l'objectif du modèle est d'évaluer si, à long terme, le système atteint un équilibre respectant les critères qui sont :
# 20% de parcelles végétalisées, dont 30% d'herbes et 70% de buissons, et garder une diversité minimale entre les deux espèces de buissons.
# # Implémentation

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
s = [160, 12, 14, 14] #200 parcelles et 50 plantées, pas d'herbes initialement parce que l'objectif final est 70% de buissons parmi la végétation, donc si on met les herbes on risque d'en avoir trop à l'équilibre
states = length(s)
patches = sum(s)

# Matrice de transitions
T = zeros(Float64, states, states)
T[1, :] = [95, 2, 1.5, 1.5] ## Barren (sol nu dominant, mais colonisation possible. Parcelle nue peut devenir herbe, buisson 1, buisson 2)
T[2, :] = [30, 50, 10, 10] ## Grass (herbes persistent, mais peuvent quand même devenir des buissons)
T[3, :] = [20, 5, 70, 5] ## Shrub1 (buisson 1 persiste)
T[4, :] = [20, 5, 5, 70] ## Shrub2 (buisson 2 persiste, mais il ne domine pas nécessairement le buisson 1, maintient de la diversité)
T

states_names = ["Barren", "Grass", "Shrub1", "Shrub2"]
states_colors = [:grey40, :orange, :teal, :purple]

# Vérification des critères

function verification_equilibre(resultat)
    final = resultat[:, end]
    Grass = final[2]
    Shrub1 = final[3]
    Shrub2 = final[4]

    vegetation = Grass + Shrub1 + Shrub2
    shrubs = Shrub1 + Shrub2

    if vegetation == 0
        return false
    end
## Critères
# Cond1: Nombre total de parcelles végétalisées
# Végétation totale environ 40 parcelles (20% de 200).
# Marge de plus ou moins 8 parcelles pour tenir compte de la variabilité stochastique.
# Cette plage reste centrée sur l'objectif de 20% et permet d'évaluer un équilibre réaliste.

    condition1 = abs(vegetation - 40) <= 8 

# Cond2: Proportion d'herbes (Grass) parmi la végétation.
# Vise 30% d'herbes (30% du 20% de 200).
# Marge de plus ou moins 0.15 pour réfléter les fluctuations stochastiques entre les simulations.
# La proportion d'herbes reste proche de la proportion cible sans exiger une valeur exacte à chaque simulation.

    condition2 = abs(Grass / vegetation - 0.3) <= 0.15

#Proportion des buissons (Shrub1 et 2) parmi la végétation.
# Vise 70% de buissons (70% du 20% de 200).
# Marge de 0.15, même logique que pour les herbes
# Proportion de buissons domine la végétation.

    condition3 = abs(shrubs / vegetation - 0.7) <= 0.15

# Diversité minimale entre les deux types de buissons
# Le buisson le moins abondant doit représenter au moins 30% du total des buissons.
# Ne peut pas tolérer de marges, on applique exactement le seuil demandé.

    condition4 = min(Shrub1, Shrub2) >= 0.30 * shrubs

    return condition1 && condition2 && condition3 && condition4
end

## #Simulations

#Nombre de simulations à effectuer
nombre_simulations = 100

#Compteur des simulations qui respectent les critères
nombre_reussites = 0

for i in 1:nombre_simulations
    resultat = simulation(T, s; stochastic=true, generations=200)
    if verification_equilibre(resultat)
        nombre_reussites += 1
    end
end

# Calcul du pourcentage de réussite (Vise au moins 80%)
pourcentage = (nombre_reussites / nombre_simulations) * 100
println("Pourcentage de réussite: ", pourcentage, "%")

## # Graphique
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
