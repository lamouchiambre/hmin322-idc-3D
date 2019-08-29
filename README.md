# HMIN322 : idc-3D

Projet *template* pour le TP sur l'insertion de données cachées dans les objets 3D.
Ce document contient la liste des activités à réaliser durant le TP.

## Pré-requis

- git
- cmake
- g++/gcc
- make
- libblas-dev liblapack-dev libxrandr-dev libxinerama-dev libxcursor-dev libxi-dev freeglut3-dev
- (sudo apt-get install <libname>)

## Préparation

- Cloner le dépôt
```sh
git clone https://github.com/TsubameDono/hmin322-idc-3D.git
```
- Compiler le projet de base. Construit automatiquement la bibliothèque **libigl** et télécharger toutes les bibliothèques externes (eigen, etc.)
```sh
$ mkdir build # create build directory
$ cd build # enter build directory
$build/ cmake .. # build libigl project (automatically load external libs)
#$ build/ cmake .. -G "Visual Studio 15 2017 Win64" # for window
$build/ make # build libraries
OR
$build/ cmake --build .
```

## Explication

**libigl** : librairie contenant des méthodes de traitement des objets 3D (IO, lissage, remaillage, filtrage, métrique)
 - igl::readXXX(.) (OFF, STL, etc.)
 - igl::writeXXX(.)
 - igl::hausdorff(.)
 - igl::is_border_vertex(.)
**libigl/opengl** : sous-librairie permettant de prototyper un outil de rendu avec igl
```cpp
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>

Eigen::MatrixXd V;
Eigen::MatrixXi F;

int main(int argc, char *argv[])
{
  // Load a mesh in OFF format
  igl::readOFF(TUTORIAL_SHARED_PATH "/bunny.off", V, F);

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.launch();
}
```
**eigen** : librairie de structures mathématiques matricielles pour la représentation
 - Eigen::MatrixXf ou Eigen::MatrixXd pour les coordonnées, Eigen::MatrixXi pour les indices
```cpp
Eigen::MatrixXd V;
Eigen::MatrixXi F;
```

## Travaux

1. Clôner le dépôt
2. Charger un maillage 3D à l'aide de la bibliothèque **libigl** en utilisant la fonction *read_triangle(.)*
3. Sauvegarder une copie du maillage 3D à l'aide de la bibliothèque **libigl** en utilisant la fonction *write_triangle(.)*
