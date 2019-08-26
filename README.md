# HMIN322 : idc-3D

Projet *template* pour le TP sur l'insertion de données cachées dans les objets 3D.
Ce document contient la liste des activités à réaliser durant le TP.

## Pré-requis

- cmake
- g++/gcc
- make

## Préparation

Permet de charger toutes les bibliothèques du projet (libigl, eigen, etc.)
```sh
git submodule update --init --recursive
```
## Travaux

1. Clôner le dépôt
2. Charger un maillage 3D à l'aide de la bibliothèque **libigl** en utilisant la fonction *read_triangle(.)*
3. Sauvegarder une copie du maillage 3D à l'aide de la bibliothèque **libigl** en utilisant la fonction *write_triangle(.)*
