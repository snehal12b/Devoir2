
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

# # Introduction

# La gestion de la végétation sous les lignes électriques à haute tension est considérée comme un défi important.
# Il faut assurer la sécurité des infrastructures et préserver aussi la biodiversité. 
# Une végétation trop dense ou composée d'arbre de grande taille peut interférer avec les lignes électriques, tandis qu'une
# absence complète de végétation peut diminuer la diversité biologique et ainsi favoriser l'érosion des sols.
# Une solution serait de créer des corridors végétalisés principalement composés d'herbes et de buissons de petite taille.
# Cela permet donc de maintenir une couverture végétale et d'éliminer les risques pour les infrastructures.
# Ce travail utilise un modèle de transition végétale pour simuler la colonisation et la succession écologique dans un corridor divisé en plusieurs parcelles.
# L'objectif est donc de déterminer une population initiale de plantes et une matrice de transition décrivant les probabilités de changement d'état des parcelles
# qui permettrait d'atteindre un équilibre écologique souhaité. En effet, le mandat indique qu'à l'équilibre, 20% des parcelles doivent être végétalisées, dont 30%
# d'herbes et 70% de buissons. Aussi, pour maintenir une diversité minimale, la variété de buisson la moins abondante doit représenter au moins 30% des parcelles 
# occupées par des buissons. À l'aide de simulations stochastiques et déterministes, nous souhaitons évaluer si les paramètres choisis permettent de respecter ces
# critères dans au moins 80% des simulations.

# # Présentation du modèle

# Le corridor étudié est représenté comme un ensemble de 200 parcelles qui peuvent se trouver dans différents états de végétation: Barren (sol nu), Grass (herbes), Shrub1 (buisson 1) et Shrub2 (buisson 2).
# Les parcelles évoluent d'un état à l'autre au fil du temps selon une matrice de transition, qui décrit les probabilités de colonisation, de persistance ou de remplacement entre les différents types de végétation.
# Le modèle utilisé est un modèle de transition markovien où l'état futur d'une parcelle dépend juste de son état actuel et des probabilités indiquées dans la matrice de transition. 
# Les transitions permettent la colonisation des parcelles nues par des herbes (Grass) ou des buissons (Shrub1 et 2), la persistance des végétaux déjà présents, et la conversion entre les différents types de végétation.
# Pour la simulation, la population initiale est composée principalement de parcelles nues (Barren), avec une petite proportion de parcelles végétalisées pour guider la succession écologique vers l'équilibre souhaité.
# En effet, pour simuler l'évolution de la végétation dans le corridor en accord avec les objectifs du mandat,nous avons choisi une population initiale composée de 160 parcelles nues (Barren), 12 parcelles d'herbes (Grass), et 14 parcelles pour chacun des deux types de buissons (Shrub1 et Shrub2).
# Donc, l'objectif du modèle est d'évaluer si, à long terme, le système atteint un équilibre respectant les critères définis dans le mandat qui sont: 20% de parcelles végétalisées, dont 30% d'herbes et 70% de buissons, tout en maintenant une diversité minimale entre les deux espèces de buissons.

# # Implémentation

# Le modèle a été créé en Julia afin de simuler l'évolution de la végétation dans un corridor composé de 200 parcelles. Les états que peuvent prendre 
# les parcelles sont les suivants: Barren (sol nu), Grass (herbes), Shrub1 (buisson 1) et Shrub2 (buisson 2). 
# Une matrice de transition a été créé pour décrire les probabilités de changement d'état d'une génération à l'autre. Une fonction d'abord vérifie que chaque ligne de la matrice 
# de transition est valide en s'assurant que la somme des probabilités de chaque ligne est égale à 1. Sinon, les valeurs sont automatiquement normalisées. 
# La population initiale est définie par un vecteur indiquant le nombre de parcelles dans chaque état. Elle est composée de 160 parcelles nues, 12 parcelles d'herbes, 
# et 14 parcelles pour chacun des deux types de buissons. La quantité d'herbe est volontairement faible pour éviter d'avoir une dominance trop importante d'herbes puisqu'on
# souhaite obtenir 70% de buissons dans la végétation.
# Il y a deux types de simulations, la simulation déterministe et la simulation stochastique. La simulation déterministe répresente la tendance moyenne du système 
# alors que la simulation stochastique permet de représenter la variabilité naturelle du système.
# Les simulations sont réalisées sur 200 générations. Pour la simulation stochastique, 100 répétitions sont réalisées pour vérifier la fiabilité des résultats et donc pour vérifier si les conditions du mandat sont respectées.
# Après chaque simulation, une fonction vérifie si quatre conditions sont respectées. La première est que le nombre total de parcelle végétalisées doit être d'environ 
# 40 parcelles. La deuxième est que les herbes doivent représenter environ 30 % de la végétation. La troisième est que 
# les buissons doivent représenter environ 70% de la végétation. La quatrième est que le type de buisson le moins abondant doit représenter au moins 30% des parcelles
# occupées par des buissons. Si les quatre conditions sont respectées alors la simulation est réussie. 
# Le pourcentage de simulations qui respectent ces conditions est ensuite calculé pour voir si le modèle atteint l'objectif de 80% de réussite.
# Finalement, un graphique est créer afin de visualiser l'évolution du nombre de parcelles au cours du temps. Les lignes pâles représentent les simulations stochastiques
# et les lignes épaisses représentent la simulation déterministe. 

# ## Packages nécessaires

import Pkg
Pkg.add("Distributions")
Pkg.add("CairoMakie")

using CairoMakie
using Distributions

import Random
Random.seed!(123456)
# # Code pour le modèle

## Vérifie que chaque ligne de la matrice de transition est valide.
## Si une ligne n'est pas valide, elle est normalisée pour que la somme soit égale à 1. Un avertissement est alors affiché pour indiquer que la ligne a été modifiée.

function check_transition_matrix!(T)
    for ligne in axes(T, 1)
        if sum(T[ligne, :]) != 1
            @warn "La somme de la ligne $(ligne) n'est pas égale à 1 et a été modifiée"
            T[ligne, :] ./= sum(T[ligne, :])
        end
    end
    return T
end

## Vérifie que les dimensions de la matrice et du vecteur états sont correctes.
## La matrice de transition doit être carrée et le nombre d'états doit correspondre à la taille de la matrice. Si ce n'est pas le cas, une erreur est levée.

function check_function_arguments(transitions, states)
    if size(transitions, 1) != size(transitions, 2)
        throw("La matrice de transition n'est pas carrée")
    end

    if size(transitions, 1) != length(states)
        throw("Le nombre d'états ne correspond psa à la matrice de transition")
    end
    return nothing
end

# Fonction simulation stochastique

function _sim_stochastic!(timeseries, transitions, generation)
    for state in axes(timeseries, 1)
        pop_change = rand(Multinomial(timeseries[state, generation], transitions[state, :]))
        timeseries[:, generation+1] .+= pop_change
    end
end

# Fonction simulation déterministe

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

# ## États, population initiale et matrice
# States
# Barren, Grass, Shrub1, Shrub2

# Population initiale
s = [160, 12, 14, 14] # 200 parcelles et 40 plantées, peu d'herbes initialement parce que l'objectif final est 70% de buissons parmi la végétation, donc si on met plus d'herbes on risque d'en avoir trop à l'équilibre.
states = length(s)
patches = sum(s)

# Matrice de transitions
T = zeros(Float64, states, states)
T[1, :] = [95, 2, 1.5, 1.5] # Barren (sol nu dominant, mais colonisation possible. Parcelle nue peut devenir herbe, buisson 1, buisson 2)
T[2, :] = [30, 50, 10, 10] # Grass (herbes persistent, mais peuvent quand même devenir des buissons)
T[3, :] = [20, 5, 70, 5] # Shrub1 (buisson 1 persiste)
T[4, :] = [20, 5, 5, 70] # Shrub2 (buisson 2 persiste, mais il ne domine pas nécessairement le buisson 1, maintient de la diversité)
T

states_names = ["Barren", "Grass", "Shrub1", "Shrub2"]
states_colors = [:grey40, :orange, :teal, :purple]

# ## Vérification des critères

function verification_equilibre(resultat)
    final = resultat[:, end] # Sélectionne la dernière colonne du résultat, qui correspond à l'état final après les générations simulées.
    Grass = final[2] # Récupère le nombre de parcelles d'herbes (Grass) à l'équilibre.
    Shrub1 = final[3] # Récupère le nombre de parcelles de buisson 1 (Shrub1) à l'équilibre.
    Shrub2 = final[4] # Récupère le nombre de parcelles de buisson 2 (Shrub2) à l'équilibre.

    vegetation = Grass + Shrub1 + Shrub2 # Total de parcelles végétalisées à l'équilibre
    shrubs = Shrub1 + Shrub2 # Total de parcelles de buissons.

    if vegetation == 0 # Vérification du cas où il n'y a aucune parcelle végétalisée.
        return false
    end

    condition1 = abs(vegetation - 40) <= 8 
    condition2 = abs(Grass / vegetation - 0.3) <= 0.15
    condition3 = abs(shrubs / vegetation - 0.7) <= 0.15
    condition4 = min(Shrub1, Shrub2) >= 0.30 * shrubs

    return condition1 && condition2 && condition3 && condition4 # On utilise &&, car && n'évalue pas les conditions suivantes si une condition est fausse, ce qui est plus efficace que & qui évalue toutes les conditions même si une est fausse.
end

## Cond1: Nombre total de parcelles végétalisées
## Végétation totale environ 40 parcelles (20% de 200).
## Marge de plus ou moins 8 parcelles pour tenir compte de la variabilité stochastique.
## Cette plage reste centrée sur l'objectif de 20% et permet d'évaluer un équilibre réaliste.

## Cond2: Proportion d'herbes (Grass) parmi la végétation.
## Vise 30% d'herbes (30% du 20% de 200).
## Marge de plus ou moins 0.15 pour réfléter les fluctuations stochastiques entre les simulations.
## La proportion d'herbes reste proche de la proportion cible sans exiger une valeur exacte à chaque simulation.

## Cond3: Proportion des buissons (Shrub1 et 2) parmi la végétation.
## Vise 70% de buissons (70% du 20% de 200).
## Marge de 0.15, même logique que pour les herbes
## Proportion de buissons domine la végétation.

## Cond4: Diversité minimale entre les deux types de buissons
## Le buisson le moins abondant doit représenter au moins 30% du total des buissons.
## Ne peut pas tolérer de marges, on applique exactement le seuil demandé pour garantir la diversité.
   
# ## Simulations réussies

# Nombre de simulations à effectuer
nombre_simulations = 100

# Compteur des simulations qui respectent les critères
nombre_reussites = 0 # Valeur initial de 0, elle augmentera à chaque fois qu'une simulation respecte les critères définis dans la fonction verification_equilibre.

for i in 1:nombre_simulations
    resultat = simulation(T, s; stochastic=true, generations=200)
    if verification_equilibre(resultat)
        global nombre_reussites += 1 # Si la simulation respecte les critères, on ajoute un au compteur de réussites. global est utilisé pour indiquer que nous faisons référence à la variable nombre_reussites définie avant la boucle.
    end
end

# Calcul du pourcentage de réussite (Vise au moins 80%)
pourcentage = (nombre_reussites / nombre_simulations) * 100
println("Pourcentage de réussite: ", pourcentage, "%")

# ## Simulations et Visualisation

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


# # Présentation des résultats

# Les résultats de la simulation montrent l'évolution du nombre de parcelles dans chacun des quatres états de végétation (Barren, Grass, Scrhub 1 et Schrub 2)   
# au cours de 200 générations. Les lignes pâles représentent les simulations stochastiques, et les lignes épaisses représentent la simulation déterministe. Un modèle déterministe 
# impose une règle fixe reliant entrées et sorties de façon prévisible. Les résultats seront toujours les mêmes pour un état intial donné. Au contraire, un  modèle stochastique 
# introduit des éléments aléatoires dans le système, ce qui permet, en programmation de simulations, de représenter des variations d’un phénomène de manière plus naturelle 
# en produisant différents résultats à chaque exécution, avec les mêmes paramètres de départ (husson2001modele).
# Au début de la simulation, la majorité des parcelles sont barren (environ 160 parcelles), alors que les parcelles végétalisées représentent environ 40 parcelles. 
# Au cours du temps, le nombre de parcelles dans chaque état évolue progressivement vers un équilibre stable. Les résultats montrent que les parcelles 
# Barren restent dominantes, tandis que les autres états occupent une proportion plus faible, mais constante. Les simulations stochastiques varient mais la tendence 
# générale reste similaire entre les simulations. La plupart convergent vers des valeurs proche de celles de la  simulation déterministe.
# Avec un taux de réussite de 80% dans les simulations stochastiques, cela indique que le modèle atteint les objectifs définis par le mandat dans la majorité des cas. 
# Ceci suggère que les paramètres choisis sont appropriés pour atteindre l'équilibre souhaité.


# # Discussion

# Les résultats des simulations montrent que le système converge vers un équilibre stable après une phase initiale de transition et après plusieurs générations. La majorités des parcelles restent Barren, tandis
# que le reste du corridor est occupé par de la végétation. À l'équilibre, le nombre de parcelles végétalisées est autour de 40 parcelles, ce qui correspond à environ 20 % des 200 parcelles.
# Cela indique que la matrice de transition choisie permet globalement de maintenir la proportion de végétation souhaitée. La répartition des herbes et buisssons demeure relativement stable
# au cours des simulations. Les herbes colonisent activement certaines parcelles nues, tandis que les buissons ont une probabilité de persistance plus élevée.
# De plus, les deux types de buissons restent présents dans des proportions relativement similaires dans la plupart des simulations, ce qui permet de maintenir la diversité minimale requise.
# La différence entre les modèles déterministes et stochastiques est importante pour interpréter les résultats. Le modèle déterministe montre la tendance moyenne de l'évolution de la végétation. Cela
# sert comme référence pour comprendre la dynamique générale. Le modèle stochastique prend en compte la variabilité des processus de colonisation des parcelles ainsi que la persistance des buissons. 
# C'est pour cela que certaines simulations peuvent montrer des dominances temporaires d'herbes ou un déséquilibre entre les deux espèces de buissons. Ces varations permetent de représenter la réalité 
# écologique, car la colonisation et succession sont influencées par des facteurs aléatoires. 
# Cependant, la nature stochastique du modèle engendre une variabilité entre les simulations. Bien que le taux de réussite soit de 80%, certains scénarios ne respectent pas exactement les critères du mandat.
# Ces échecs ont lieu principalement lorsque qu'il y a une dominance d'herbes (Grass) ou un manque d'équilibre entre les deux types de buissons (Shrub1 et Shrub2).
# Par exemple, dans certaines simulations, les herbes peuvent coloniser plus rapidement que les buissons, ce qui peut entraîner une proportion d'herbes plus élevée que les 30% souhaités.
# Les résultats, montrent la sensibilité du système. En effet, même si la matrice de transition est conçue pour favoriser un équilibre respectant les critères du mandat, la variabilité naturelle du système
# peut mener à des résultats différents d'une simulation à l'autre.
# Ce modèle permet donc de combiner la prédicitibilité (déterminisme) et réalisme écologique (stochastique) dans la modélisation de la végétation sous lignes à haute tension.

 # # Bibliographie
