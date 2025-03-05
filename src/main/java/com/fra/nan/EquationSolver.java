package com.fra.nan;

import java.util.Scanner;

// Implémentation de l'interface Imat
public class EquationSolver implements Imat {
    @Override
    public double[] resoudreEquation(double a, double b, double c) {
        if (a == 0) {
            double x = - c / b; 
            System.out.println("Solution unique : " + x);
            return new double[]{x};
            // throw new IllegalArgumentException("Ce n'est pas une équation du second degré.");
        }
        
        double delta = b * b - 4 * a * c;
        
        if (delta > 0) {
            double x1 = (-b + Math.sqrt(delta)) / (2 * a);
            double x2 = (-b - Math.sqrt(delta)) / (2 * a);
            System.out.println("Solutions : " + x1 + " ... " + x2);
            return new double[]{x1, x2};
        } else if (delta == 0) {
            double x = -b / (2 * a);
            System.out.println("Solution double : " + x);
            return new double[]{x};
        } else {
            System.out.println("Pas de solution");
            return new double[]{}; // Pas de solution réelle
        }
    }
    
    
    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);
        
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
        
        scanner.close();
    }
}

/* class Main {
    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);
        
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
        
        scanner.close();
    }
}
*/