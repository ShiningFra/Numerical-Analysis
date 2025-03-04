/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package com.fra.nan;

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
        0.1,        // 3 - Petit réel positif pour a
        -0.1,       // 4 - Petit réel négatif pour a
        1000.0,     // 5 - Grand réel positif pour a
        -1000.0,    // 6 - Grand réel négatif pour a
        0.0,        // 7 - Cas limite pour b
        2.0,        // 8 - Cas normal positif pour b
        -2.0,       // 9 - Cas normal négatif pour b
        0.5,        // 10 - Petit réel positif pour b
        -0.5,       // 11 - Petit réel négatif pour b
        100.0,      // 12 - Grand réel positif pour b
        -100.0,     // 13 - Grand réel négatif pour b
        0.0,        // 14 - Cas limite pour c
        1.0,        // 15 - Cas normal positif pour c
        -1.0,       // 16 - Cas normal négatif pour c
        0.1,        // 17 - Petit réel positif pour c
        -0.1,       // 18 - Petit réel négatif pour c
        1000.0,     // 19 - Grand réel positif pour c
        -1000.0     // 20 - Grand réel négatif pour c
    };

    // Fonction pour associer les lettres de l'alphabet aux entrées du tableau
    public static char getEntryLetter(int index) {
        return (char) ('A' + index);
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
