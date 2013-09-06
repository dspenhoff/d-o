import java.io.*;
import java.util.List;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Random;

public class Solver {
  
  int N;
  Point[] points; 
  int[] tour;
  double tourLength;     // the objective function
  boolean[] visited;      // has a point been visited yet on the tour? 
  private double temp, minTemp;
  private double cooling;
  private int tempCount, tempCountMax; 
  private Random rng;
    
  /**
  * The main class
  */
  public static void main(String[] args) {
    try {
      Solver solver = new Solver(args);
    } catch (IOException e) {
        e.printStackTrace();
    }
  }
    
  /**
  * Read the instance, solve it, and print the solution in the standard output
  */
  public Solver(String[] args) throws IOException {
    String fileName = null;
        
    // get the temp file name
    for(String arg : args){
      if(arg.startsWith("-file=")){
        fileName = arg.substring(6);
      } 
    }
    if(fileName == null)
      return;
        
    // read the lines out of the file
    List<String> lines = new ArrayList<String>();

    BufferedReader input =  new BufferedReader(new FileReader(fileName));
    try {
      String line = null;
      while (( line = input.readLine()) != null){
        lines.add(line);
      }
    }
    finally {
      input.close();
    }
              
    // parse the data in the file
    String[] firstLine = lines.get(0).split("\\s+");
    N = Integer.parseInt(firstLine[0]);
    points = new Point[N];

    for(int i = 0; i < N; i++){
      String line = lines.get(i+1);
      String[] parts = line.split("\\s+");
      points[i] = new Point(Double.parseDouble(parts[0]), Double.parseDouble(parts[1]), i);
    } 
    

    rng = new Random();
    double minLength = Double.MAX_VALUE;
    int minStart = 0;
    for (int i = 0; i < N; i++) {
      solver(rng.nextInt(N));
      if (tourLength < minLength) {
        minLength = tourLength;
        minStart = i;
        //printSolution(false);
      }
            
    }  
    
    // display the final solution
    solver(minStart);
    if (!feasibleTour()) System.out.println("infeasible");
    printSolution(true);
    
  }
  
  private void solver(int start) {
    nnGreedy(start);

    boolean refine = true;
    double prevTourLength = tourLength;
    while(refine) {
      // set up the simulated annealing
      temp = .9999;
      cooling = .99;
      minTemp = .000001; 
      tempCountMax = 100;
      tempCount = 0;
      twoOpt();
      if (tourLength < prevTourLength) {
        prevTourLength = tourLength;
      } else {
        refine = false;
      }
    }


    prevTourLength = tourLength;
    refine = true;
    while(refine) {
      threeOpt();
      if (tourLength < prevTourLength) {
        prevTourLength = tourLength;
      } else {
        refine = false;
      }
    }


    //twoOpt(); threeOpt(); printSolution(false);
    //twoOpt(); threeOpt(); printSolution(false);
    //twoOpt(); threeOpt(); printSolution(false);
    //twoOpt(); threeOpt(); printSolution(false);
    //twoOpt(); threeOpt(); printSolution(false);

    //threeOpt();
    //threeOpt();
    
    refine = true;
    prevTourLength = tourLength;
    while(refine) {
      twoOpt();
      if (tourLength < prevTourLength) {
        prevTourLength = tourLength;
      } else {
        refine = false;
      }
    }

  }
  
  private double dist(Point p, Point q) {
    double tx = p.x() - q.x();
    double ty = p.y() - q.y();
    return Math.sqrt(tx * tx + ty * ty);
  }
  
  private void nnGreedy(int first) {
    // greedy solution: visit next closest point (nearest neighbor)
    tour = new int[N];
    boolean[] visited = new boolean[N];
    tour[0] = first;
    visited[first] = true;
    int prev = first;
    tourLength = 0;
    for (int i = 1; i < tour.length; i++) {
      int next = -1;
      double mind = Double.MAX_VALUE;
      for (int j = 0; j < N; j++) {
        if (!visited[j]) {
        
          // is point j outside the outer enclosing box for prev, ie not closer than min point?
          if (points[j].x() > points[prev].x() + mind || points[j].x() < points[prev].x() - mind) continue;
          if (points[j].y() > points[prev].y() + mind || points[j].y() < points[prev].y() - mind) continue;
                    
          double d = dist(points[prev], points[j]);
          if (d < mind) {
            mind = d;
            next = j;
            
          }          
        }
      }
      tour[i] = next;
      visited[next] = true;
      tourLength += mind;
      prev = next;      
    }    

    tourLength += dist(points[tour[N-1]], points[tour[0]]);    
  }
  
  private void twoOpt() {
    // 2-opt local improvement
       
    for (int i = 0; i < N-1; i++) {      
      for (int j = i+1; j < N; j++) {
         
        // can the tour be improved by swapping <i, i+1>, <j, j+1>?
        int ip1 = ((i + 1) < N) ? i + 1 : 0;
        int jp1 = ((j + 1) < N) ? j + 1 : 0;
        double di = dist(points[tour[i]], points[tour[ip1]]);   // length from i
        double dj = dist(points[tour[j]], points[tour[jp1]]);
        double newdi = dist(points[tour[i]], points[tour[j]]);  // new length from i
        double newdj = dist(points[tour[ip1]], points[tour[jp1]]);
        double delta = - (di + dj) + newdi + newdj;

        if (delta < 0) {
        //if (delta < 0 || anneal(delta)) {
          
          int[] temp = new int[j-i];
          for (int k = 0; k < j-i; k++) temp[k] = tour[j-k];
          for (int k = 0; k < j-i; k++) tour[i+1+k] = temp[k];
          tourLength += delta;
        }
      }
    }    
  }
  
  private boolean anneal(double delta) {
    // cool the temp?
    if (++tempCount > tempCountMax) {
      temp *= cooling;
      if (temp < minTemp) temp = minTemp;
      tempCount = 0;
    }
    return (rng.nextDouble() < Math.exp(- delta / temp));
  }
  
  private void threeOpt() {
    // 3-opt local improvement
    for (int i = 0; i < N-3; i++) {      
      for (int j = i+1; j < N-2; j++) {
        for (int k = j+1; k < N-1; k++) {
           
           // can the tour be improved by swappping <i, i+1>, <j, j+1>, <k, k+1>?
           int ip1 = ((i + 1) < N) ? i + 1 : 0;
           int jp1 = ((j + 1) < N) ? j + 1 : 0;
           int kp1 = ((k + 1) < N) ? k + 1 : 0;
           double di = dist(points[tour[i]], points[tour[ip1]]);   // current segment lengths
           double dj = dist(points[tour[j]], points[tour[jp1]]); 
           double dk = dist(points[tour[k]], points[tour[kp1]]);          
           double newdi1 = dist(points[tour[i]], points[tour[jp1]]);  // new segment lengths
           double newdj1 = dist(points[tour[j]], points[tour[kp1]]);                    
           double newdk1 = dist(points[tour[k]], points[tour[ip1]]);
           double delta1 = - (di + dj + dk) + newdi1 + newdj1 + newdk1;
           double newdi2 = dist(points[tour[i]], points[tour[jp1]]);  // alternate new segment lengths
           double newdj2 = dist(points[tour[k]], points[tour[j]]);                    
           double newdk2 = dist(points[tour[kp1]], points[tour[ip1]]);
           double delta2 = - (di + dj + dk) + newdi2 + newdj2 + newdk2; 

           if (delta1 < 0 || delta2 < 0) {
             //System.out.println("delta1 = " + delta1);
             if (delta1 < delta2) {
               int[] tempi = new int[j-i];
               int[] tempj = new int[k-j];
               for (int l = 0; l < k-j; l++) tempj[l] = tour[j+1+l];
               for (int l = 0; l < j-i; l++) tempi[l] = tour[i+1+l];
               for (int l = 0; l < k-j; l++) tour[i+1+l] = tempj[l];
               for (int l = 0; l < j-i; l++) tour[i+k-j+1+l] = tempi[l];             

               //System.out.println("\n prev length, delta = " + tourLength + ", " + delta1);
               tourLength += delta1;
               //System.out.println(tourLength);
               //System.out.println(feasibleTour());
               // printSolution(false);
               //return;               
             } else {
                int[] tempi = new int[j-i];
                int[] tempj = new int[k-j];
                for (int l = 0; l < k-j; l++) tempj[l] = tour[j+1+l];
                for (int l = 0; l < j-i; l++) tempi[l] = tour[j-l];
                for (int l = 0; l < k-j; l++) tour[i+1+l] = tempj[l];
                for (int l = 0; l < j-i; l++) tour[i+k-j+1+l] = tempi[l];  
                tourLength += delta2;             
             }
             //System.out.println("i, i, k = " + i + ", " + j + ", " + k);
             //System.out.println(tourLength);
             //System.out.println(feasibleTour());
             //printSolution(true);
           } 
        } 
      }
    }   
  }
  
  private boolean feasibleTour() {
    // is the tour feasible, eg is each point in the tour exactly one
    // sum of tour points should equal sum of N-1 integers
    int count = 0;
    for (int i = 0; i < N; i++) count += tour[i];
    return (count == ((N) * (N-1) / 2));
  }
  
  
  private double printSolution(boolean printTour) {
    // print the solution in the specified output format        
    double length = 0;
    String s = tour[0] + " ";
    for (int i = 0; i < tour.length-1; i++) {
      length += dist(points[tour[i]], points[tour[i+1]]);
      s += tour[i+1] + " ";
    }
    length += dist(points[tour[N-1]], points[tour[0]]);
    System.out.println(length + " 0");
    if (printTour) System.out.println(s);
    return length;
  }
  
  private double length() {
    double l = 0;
    for (int i = 0; i < tour.length-1; i++) {
      l += dist(points[tour[i]], points[tour[i+1]]);
    }
    return l + dist(points[tour[N-1]], points[tour[0]]);
    
  }

  private class Point {
    private double x, y;
    int point;
    public Point (double a, double b, int i) {
      x = a;
      y = b;
      point = i;
    }
    public double x() { return x; }
    public double y() { return y; }
    public int point() { return point; }
    public String toString() {
      return "(" + x + ", " + y + ")";
    }
  }
  
  private class TourSegment implements Comparable<TourSegment> {
    private int from, to;   // indices into the points array
    private int index;      // from's index in tour
    private double distance;
    private TourSegment(int p, int q, int i) {
      from = p;
      to = q;
      index = i;
      distance = dist(points[from], points[to]);
    }
    public int compareTo(TourSegment that) {
      if (this.distance > that.distance) return 1;
      if (this.distance < that.distance) return -1;
      return 0;
    }
    private int from() { return from; }
    private int to() { return to; }
    private void setTo(int q) {
      to = q;
      distance = dist(points[from], points[to]);
    }
    private double distance() { return distance; }
    public String toString() {
      return "segment from: " + from + " to: " + to + ", distance: " + distance;
    }

    
  }
  
}