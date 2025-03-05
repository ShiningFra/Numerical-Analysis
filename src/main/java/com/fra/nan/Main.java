package com.fra.nan;

import java.util.Scanner;

public class Main {
    public static void main(String[] args) {
        /*Scanner scanner = new Scanner(System.in);
        
        System.out.print("Entrez le coefficient a : ");
        double a = scanner.nextDouble();
        System.out.print("Entrez le coefficient b : ");
        double b = scanner.nextDouble();
        System.out.print("Entrez le coefficient c : ");
        double c = scanner.nextDouble();
        
        Imat solver = new EquationSolver();
        double[] solutions = solver.resoudreEquation(a, b, c);
        
        if (solutions.length == 2) {
            System.out.println("Les solutions sont : x1 = " + solutions[0] + ", x2 = " + solutions[1]);
        } else if (solutions.length == 1) {
            System.out.println("La solution unique est : x = " + solutions[0]);
        } else {
            System.out.println("Pas de solution réelle.");
        }
        
        scanner.close();*/
        EquationTest.testResoudreEquation();
    }
}
