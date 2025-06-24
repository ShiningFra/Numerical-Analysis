# Analyse Numérique

Ce projet contient des solveurs pour l'équation différentielle `-u'' = f` (en 1D) et `-Δu = f` (en 2D) en utilisant les méthodes des Différences Finies et des Volumes Finis.

## Structure du Projet

Le code source Java se trouve dans `src/main/java/`. Chaque méthode a son propre package :
*   `FiniteDifference1` : Différences Finies 1D
*   `FiniteVolume1` : Volumes Finis 1D
*   `FiniteDifference2` : Différences Finies 2D
*   `FiniteVolume2` : Volumes Finis 2D

Chaque package contient une classe principale avec une méthode `main` qui lance une interface graphique Swing pour cette méthode spécifique (par exemple, `FiniteDifference1.ODEFiniteDifference`).

## Prérequis

*   JDK (Java Development Kit) installé (par exemple, version 8 ou ultérieure).
*   Accès à un terminal ou une ligne de commande.

## Compilation

1.  Ouvrez un terminal ou une ligne de commande.
2.  Naviguez jusqu'à la racine du projet (où se trouve le dossier `src`).
3.  Compilez tous les fichiers Java.

- Créez un répertoire pour les classes compilées (par exemple, `out`):
    ```bash
    mkdir out
    ```
- Compilez les fichiers `.java` pour les DF et VF, 1D et 2D en spécifiant le répertoire de sortie et l'encodage UTF-8 pour une gestion correcte des accents. Depuis la racine du projet :
    ```bash
    javac -encoding UTF-8 -d out src/main/java/FiniteDifference1/ODEFiniteDifference.java src/main/java/FiniteVolume1/ODEFiniteVolume.java src/main/java/FiniteDifference2/ODEFiniteDifference.java src/main/java/FiniteVolume2/ODEFiniteVolume.java
   ```

- Compilez les fichiers `.java` pour le TP sur Jenny en spécifiant le répertoire de sortie et . Depuis la racine du projet :
    ```bash
    javac -encoding UTF-8 -d out src/main/java/Jenny/Main.java src/main/java/Jenny/EquationSolver.java src/main/java/Jenny/EquationTest.java src/main/java/Jenny/Constants.java src/main/java/Jenny/Imat.java
    ```
  (Assurez-vous que le fichier test_cases.txt est bien à la racine du projet ou ajustez le chemin dans EquationTest.java si nécessaire. Actuellement, il le cherche à la racine : new FileReader("test_cases.txt")).


## Exécution

Après la compilation, vous pouvez exécuter chaque programme en spécifiant le classpath et la classe principale.

**Si compilé manuellement dans `out` (depuis la racine du projet) :**

*   **Jenny:**
    ```bash
    java -cp out Jenny.Main
    ``` 
    jenny.exe est un fichier exécutable Windows. Pour l'exécuter :

    - Sous Windows :
        - Ouvrez une invite de commandes (cmd) ou PowerShell.
        - Naviguez jusqu'au répertoire src/main/java/Jenny/.
            ```bash
                cd src\main\java\Jenny
            ```
        - Exécutez-le directement : **jenny.exe**
    - Sous Linux/macOS :
        - Les fichiers .exe ne s'exécutent pas nativement. Vous pourriez avoir besoin d'un outil comme Wine pour tenter de l'exécuter, mais sa compatibilité n'est pas garantie.

*   **Volumes Finis 1D:**
    ```bash
    java -cp out FiniteVolume1.ODEFiniteVolume
    ```
*   **Différences Finies 2D:**
    ```bash
    java -cp out FiniteDifference2.ODEFiniteDifference
    ```
*   **Volumes Finis 2D:**
    ```bash
    java -cp out FiniteVolume2.ODEFiniteVolume
    ```


## Utilisation des Interfaces Graphiques

*   **Programmes 1D (`ODEFiniteDifference` dans `FiniteDifference1`, `ODEFiniteVolume` dans `FiniteVolume1`) :**
    *   Une fenêtre s'ouvrira.
    *   **Sélectionnez la "Solution exacte u(x)" de référence** (pour laquelle `f(x) = -u''(x)` sera calculé). Les options sont `u(x) = sin(πx)` et `u(x) = x³`.
    *   **Sélectionnez le "Type d'équation"**. Pour cette tâche, utilisez toujours `TYPE3 ("-u'' = f")`.
    *   Les conditions aux limites sont fixées dans le code à `u(0)=0` et `u(1)=1`.
    *   Cliquez sur **"Résoudre et Afficher Solution"** pour voir des graphiques de la solution numérique vs. analytique pour N=10, 20, 40, 80, 160, 320. Chaque graphique s'ouvre dans un onglet. La console affichera l'erreur L∞.
    *   Cliquez sur **"Afficher Courbe d'Erreur"** pour calculer les solutions pour N=10, 20, 40, 80, 160, 320. La console affichera l'erreur L∞ et l'ordre de convergence. Une nouvelle fenêtre s'ouvrira avec un graphique log-log de l'erreur L∞ en fonction de N.

*   **Programmes 2D (`ODEFiniteDifference` dans `FiniteDifference2`, `ODEFiniteVolume` dans `FiniteVolume2`) :**
    *   Une fenêtre s'ouvrira.
    *   **Sélectionnez la "Sol. exacte u(x,y)" de référence** (pour laquelle `f(x,y) = -Δu(x,y)` sera calculé). Options : `u(x,y) = sin(πx)sin(πy)` et `u(x,y) = x³y³`.
    *   Les conditions aux limites de Dirichlet sont prises à partir de la solution exacte choisie sur les bords du domaine `(0,1)x(0,1)`.
    *   Cliquez sur **"Résoudre et Afficher"**.
    *   La zone de texte affichera l'erreur L∞ et l'ordre de convergence pour N x N avec N = 10, 20, 40, 80, 160, 320.
    *   La partie inférieure affichera trois heatmaps pour N=80 : solution numérique, analytique, et erreur absolue.

## Interprétation des Résultats

*   **Erreur L∞ :** Différence maximale absolue entre solution numérique et analytique. Plus petite = mieux.
*   **Ordre de Convergence Numérique :** Théoriquement 2 pour les schémas utilisés. Si N double, l'erreur L∞ devrait être divisée par 4.
*   **Graphiques (1D) :**
    *   Solution : Numérique (rouge) vs Analytique (bleu). Devraient se superposer pour N grand.
    *   Courbe d'erreur (log-log) : Droite de pente ~-2 si ordre 2.
*   **Heatmaps (2D) :**
    *   Numérique et Analytique : Devraient se ressembler. Couleurs : Bleu (min) -> Rouge (max).
    *   Erreur Absolue : Devrait être majoritairement bleue (valeurs faibles).


    
