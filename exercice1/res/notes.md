# Links

[Slides](https://moodle.epfl.ch/pluginfile.php/77921/mod_resource/content/266/Presentations/PhysNum24Week01Intro.pdf)

[Cours][https://moodle.epfl.ch/pluginfile.php/77901/mod_resource/content/16/physnumbook_24_2.pdf]

# Cours 1

## Discretisation

Monde réel:
- continu $\frac{dy}{dt} = f(y,t)$
- solution exacte
- calcul différentiel et intégral

Monde numérique:
- discret $\frac{y_{n+1}-y_n}{t_{n+1}-t_n} \approx f(y_n,t_n)$
- solutions approchées
- nombres finis
- opérations arithmétiques

Important: vérifier & valider
- Convergence
- Stabilité
- Propriétés physiques (comparer avec solution exacte si possible)

## Intégration, différentiation

Dérivée: $$\frac{dy}{dt} = \underbrace{\frac{y(t+\Delta t) - y(t)}{\Delta t}}_{\textrm{Erreur d'arrondi}} + \underbrace{O(\Delta t)}_{\textrm{Erreur de troncature}}$$

## Différences finies et devs limités

3 types pour trouver les dérivées:
1. Diff finie progressive
2. Diff finie centrée
3. ... rétrograde

En utilisant les différences centrées on a des erreurs d'ordre plus grand -> plus efficace

# Evolution temporelle

Schémas numériques: ordre de convergence: pente dans un loglog, dans la limite des petits $\Delta t$

Constats:
- $\Delta t$ trop grand: ca oscille: non physique
- $\Delta t$ encore plus grand: ca oscille et explose: instable
- $\Delta t$ petit: convergence

Euler implicite: current iteration ($y_{n+1}$). Solve for $y_{n+1}$, puis résoudre par itérations (point fixe) 

Euler explicite: previous iteration ($y_n$)

Euler semi-implicite: faire la moyenne des 2 mdr