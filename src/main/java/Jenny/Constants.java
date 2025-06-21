/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package Jenny;

/**
 *
 * @author Roddier
 */
public class Constants {
    // Tableau simple avec les constantes
    public static final double[] CONSTANTS = {
        0.0,        // 0 - Cas limite pour a
        1.0,        // 1 - Cas normal positif pour a
        -1.0,       // 2 - Cas normal négatif pour a
        0.001,        // 3 - Petit réel positif pour a
        -0.001,       // 4 - Petit réel négatif pour a
        1000.0,     // 5 - Grand réel positif pour a
        -1000.0     //6 - Grand réel négatif pour a
    };

    // Fonction pour associer les lettres de l'alphabet aux entrées du tableau
    public static char getEntryLetter(int index) {
        return (char) ('a' + index);
    }
    
    public static double getValue(char letter) {
        return CONSTANTS[((int)letter) - ((int)'a')];
    }

    public static void main(String[] args) {
        // Affichage des constantes avec les lettres associées
        for (int i = 0; i < CONSTANTS.length; i++) {
            System.out.printf("Lettre : %s, Valeur : %s%n", getEntryLetter(i), CONSTANTS[i]);
        }
    }
}
