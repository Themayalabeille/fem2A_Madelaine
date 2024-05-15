# Éléments Finis en C++

#### Objectif

Écrire [un code aux éléments finis en C++](course/code.md) pour résoudre des 
[problèmes de Poisson](course/poisson.md) 2D avec des éléments triangulaires
linéaires. 

#### Contact

Paul Cupillard : paul.cupillard@univ-lorraine.fr

#### Liens utiles

Cours sur [Arche](http://arche.univ-lorraine.fr/course/view.php?id=61482)

Vidéo de Gilbert Strang : [Finite element method](https://www.youtube.com/watch?v=WwgrAH-IMOk)

Cours de Grégoire Allaire : [Approximation numérique et optimisation](http://www.cmap.polytechnique.fr/~allaire/map411/polycopie-map411.pdf)

Générateurs de maillages triangulaires : [Gmsh](http://gmsh.info/),[Triangle](https://www.cs.cmu.edu/~quake/triangle.html)

# Eléments d'utilisation du code

#### Tests

Les tests executés par le programme peuvent être choisi dans la fonction run_tests dans le fichier main.cpp. Les tests true sont executés et les test false ne le sont pas.

- ***t_Quadrature*** : teste les méthodes de la classe *quadrature*
- ***t_elem_mapping*** : teste les méthodes de la classe *ElementMapping* (sur le bord 4 avec le point (0.2, 0) et sur le triangle 4 avec le point (0.2, 0.4))
- ***t_shape_function*** : teste les méthodes de la classe *ShapeFunction*
- ***t_assemble_Ke*** : teste l'assemblage des différentes matrices et vecteurs dans le maillage square.mesh


#### Simulations

La simulation à executer est selectionnable dans la fonction run_simus du fichier main.cpp. Dans le fichier main.cpp, dans la fonction run_simus, est selectionnable la grille sur laquelle est executé la simulation en changeant le fichier .mesh ainsi que le degré de quadrature et le terme source pour la simulation avec condition de Dirichlet pur. Par défault, le degré de quadrature est demandé en console lors de l'execution d'une simulation.
Le nom du fichier solution sauvegardé peut être changé dans la fonction correspondant à la simulation choisi (dans le fichier simu.h) De même, il est possible de choisir l'option de calcul de l'erreur avec la solution analytique pour le sinus_bump.

- ***simu_pure_dirichlet*** : Simulation du problème de Dirichlet pur
- ***simu_pure_dirichlet_source*** : Simulation du problème de Dirichlet avec un terme source unitaire ou de type sinus_bump
- ***simu_neumann*** : Simulation du problème de Dirichlet et Neumann
- ***simu_mug*** : Simulation du problème du mug sur *mug_0_5.mesh* 
