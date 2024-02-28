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

# Cours 2

## Ordre de convergence

Ordre de convergence, n, tq
$$y_{n+m}(t_\text{fin}) = \underbrace{y_\text{exact}(t_\text{fin})}_{\text{valeur convergée}} + c_n(\Delta t)^n + c_{n+1}(\Delta t)^{n+1} + \cdots$$
$$c_m = 0, \, \forall m < n$$

Doit tendre vers 0 pour $\Delta t \rightarrow 0$

$n_\text{steps} \propto \frac{1}{\Delta t}$

$ \frac{dy}{dt} = f(y, t) $ (1)

## Analyse de stabilité (Von Neumann)

Euler explicite: $y_{n+1} = y_n + f(y_n,t_n) \Delta t$ (2)

Def erreur $e_n$:  $y_n = y(t_n) + e_n$ (3)

On combine (2) et (3)

$$ \underbrace{y(t_{n+1})}_{y_n(t_n)+\underbrace{\frac{dy}{dt}(t_n)}_{f(y(t_n),t_n)} \Delta t + O((\Delta t)^2)} + e_{n+1} = y(t_n)+e_n + \underbrace{f(y(t_n)+e_n,t_n)}_{f(y(t_n),t_n)+\frac{\partial f}{\partial y}(y(t_n),t_n)e_n + O(e_n^2)}$$

$$ \implies e_{n+1}=\underbrace{(1+\frac{\partial f}{\partial y}(t(t_n), t_n)\Delta t)}_{\text{Gain de l'erreur, }G} e_n $$

Si $|G| > 1$: instable, si $|G|\le 1$: stable

### Cas de la désintégration

$$f(y,t)=-\gamma y$$

$$\frac{\partial f}{\partial y} = -\gamma |G| = |1 - \gamma \Delta t|$$

$\rightarrow$ Stable su $0 < \gamma \Delta t < 2$, instable sinon.

### Généralisation

$ \frac{d\underline{y}}{dt} = \underline{f}(\underline{y}, t) $

Euler explicite: $\underline{y}_{n+1} = \underline{y}_n + \underline{f}(\underline{y}_n,t_n) \Delta t$

Def erreur $\underline{e}_n$:  $\underline{y}_n = \underline{y}(t_n) + \underline{e}_n$

$$ \underline{e}_{n+1}=\underbrace{(\mathbb{I}+\frac{\partial \underline{f}}{\partial \underline{y}}(\underline{y}(t_n), t_n)\Delta t)}_{\text{Matruce Gain de l'erreur, } \underline{\underline{G}}} \underline{e}_n $$

Add vectors, vectors everywhere...

### Cas oscillateur harmo

$$ \frac{d}{dt} \underbrace{\left(\begin{matrix} x\\ v \end{matrix}\right)}_{\underline{y}} = 
\underbrace{\underbrace{\left(\begin{matrix} 0 & 1 \\ -\frac{k}{m} & 0 \end{matrix}\right)}_{\underline{\underline{M}}}
\underbrace{\left(\begin{matrix} x\\ v \end{matrix}\right)}_{\underline{y}}}_{\underline{f}(\underline{y})} $$

On cherche les valeurs propres
$$\lambda = 1 \pm i \sqrt{\frac{k}{m}} \Delta t$$
$$\implies |\lambda| > 1, \, \forall \Delta t$$

De l'analyse de la stabilité du Von Neumann:
Si $|G| > 1$: instable, si $|G|\le 1$: stable $\iff$ $\exists i \, \text{t.q. } |\lambda_i| > 1$: instable, si $|\lambda_i| \le 1| \,\forall i$

Toujours instable avec Euler explicite!

## Solution analytique approximative des équations discrètes (pour l'oscillateur harmonique)

- Euler: $\underline{y}_{n+1} = (\mathbb{I} + \underline{\underline{M}} \Delta t) \underline{y}_n$
- Cherche solution: $\underline{y}_n = \underline{\underline{A}}_{\in \mathbb{C}^2} e^{i \tilde{\omega} t_n}$, où $\tilde{\omega} = \omega + i\gamma$ $$\underline{\underline{A}} e^{i \tilde{\omega} t_n} e^{i \tilde{\omega} \Delta t} = (\mathbb{I} + \underline{\underline{M}} \Delta t) \underline{\underline{A}} e^{i \tilde{\omega} t_n}$$
On fait l'hypothèse $|\tilde{\omega}| \Delta t \ll 1$ petit paramètre d'ordre $\sim \epsilon$ (si ...) $$\implies e^{i \tilde{\omega} \Delta t} = 1 + i \tilde{\omega} \Delta t - \frac{1}{2} \tilde{\omega}^2 (\Delta t)^2 + O(\Delta t)^3$$ $$\implies [(1 + i \tilde{\omega} \Delta t - \frac{1}{2} \tilde{\omega}^2 (\Delta t)^2) \mathbb{I} - \mathbb{I} - \underline{\underline{M}} \Delta t] \underline{\underline{A}} = 0$$ $$ \implies \det[\dots] = 0$$

Après calcul du déterminant et méthode de résolution "Ordre par ordre"

$$ \tilde{\omega} = \tilde{\omega}^{(0)} + \tilde{\omega}^{(1)} + \tilde{\omega}^{(2)} + \cdots$$
- Ordre 0: $(\tilde{\omega}^{(0)})^2 = \frac{k}{m} \implies \tilde{\omega}^{(0)} = \pm \sqrt{\frac{k}{m}}$
- Ordre 1: $ (\tilde{\omega}^{(0)} + \tilde{\omega}^{(1)})^2(1+ i\tilde{\omega}^{(0)} \Delta t + i \tilde{\omega}^{(1)} \Delta t + \cdots) = \frac{k}{m}$ $$\iff ({\tilde{\omega}^{(0)}}^2 + 2 \tilde{\omega}^{(0)} \tilde{\omega}^{(1)} + {\tilde{\omega}^{(1)}}^2)(1+ i\tilde{\omega}^{(0)} \Delta t + i \tilde{\omega}^{(1)} \Delta t + \cdots) = \frac{k}{m}$$ $$\implies \tilde{\omega}^{(1)} = i \tilde{\omega}^{(0)} \frac{\Delta t}{2} = -i\frac{k}{m} \frac{\Delta t}{2}$$

Solution finale du style $y_n = e^{\frac{k \Delta t}{m 2} t_n} |A| \cos(\sqrt{\frac{k}{m}} t_n + \varphi)$

On trouve donc un taux de croissance prop à $\Delta t$