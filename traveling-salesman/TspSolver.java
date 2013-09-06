
import java.util.Random;

public class TspSolver {
  

  /*
  * solves traveling salesman problem by creating a feasible and optimal* tour
  * that visits each point (x, y) exactly once and returns to the starting point
  *
  * if fixStartingPoint is true, point (x[0], y[0]) will be the fixed starting point
  * unless randomStart is specified
  */
  
  private int N;
  private Point[] points;
  private boolean fixed;
  private boolean random;
  private boolean three;
  private boolean anneal;
  private int[] tour;
  private int[] map;
  private double tourLength;
  private boolean[] visited;
  private Random rng;
  private Stopwatch timer;
  private double time;
  
  
  // parameters for simulated annealing
  double temperatureStart = 1000000;
  double temperature = temperatureStart;
  double temperatureMin = 0.0001;
  double temperatureDelta = 1 - temperatureMin;
  double temperatureCount = 0;
  double temperatureChange = 1000;

  public TspSolver(double[] x, double[] y, boolean randomStart, boolean fixedStart, boolean doThreeOpt, boolean doSimulatedAnnealing)  {
    
    timer = new Stopwatch();
    
    // setup points
    N = x.length;
    points = new Point[N];
    for (int i = 0; i < N; i++) points[i] = new Point(x[i], y[i], i);
    
    fixed = fixedStart;
    random = randomStart;
    if (random) fixed = true;
    three = doThreeOpt;
    anneal = doSimulatedAnnealing;
    
    rng = new Random();
    
    System.out.println(timer.elapsedTime());
    
    int minStart = 0;
    
    fixed = true;
    random = false;
    anneal = false;
    three = true;
    minStart = 1765;
    
    double minLength = Double.MAX_VALUE;
    if (!fixed) {
      for (int i = 0; i < N; i++) {
        nnGreedy(i);
        System.out.println("start, time = " + i + ", " + timer.elapsedTime());
        if (tourLength < minLength) {
          minLength = tourLength;
          minStart = i;
          System.out.println(">>>>> min start, length = " + minStart + ", " + minLength);
          
        }          
      }      
    } 
    
    if (random) minStart = rng.nextInt(N);

    // create the (final) solution from the best greedy starting point
    solver(minStart);
    
    time = timer.elapsedTime();

  }
  

  private void solver(int start) {
    nnGreedy(start);
    System.out.println(length());
    System.out.println(timer.elapsedTime());
    
    boolean refine = true;
    double prevTourLength = tourLength;
    
    while(refine) {
      System.out.println("refining");
      twoOpt();
      System.out.println(length());
      System.out.println(timer.elapsedTime());
      
      randomTwoOpt(100000);
      System.out.println(length());
      System.out.println(timer.elapsedTime());
      
      if (tourLength < prevTourLength) {
        prevTourLength = tourLength;
      } else {
        refine = false;
      }
    }

    if (three) {
      System.out.println("three");
      System.out.println(timer.elapsedTime());
      
      randomThreeOpt(100000);
      System.out.println(length());
      System.out.println(timer.elapsedTime());

      refine = true;
      while (refine) {
        threeOpt();
        System.out.println(length());
        System.out.println(timer.elapsedTime());
        if (tourLength < prevTourLength) {
          prevTourLength = tourLength;
        } else {
          refine = false;
        } 
      } 
       
    }    
  }
  
  public TspSolver() {}     // null constructor

  private void nnGreedy(int first) {
    // greedy solution: visit next closest point (nearest neighbor)
    tour = new int[N];
    visited = new boolean[N];
    add(0, first);
    for (int j = 1; j < N; j++) add(j, closestUnvisited(tour[j-1]));
    tourLength = length();    
  }
  
  private void add(int i, int j) {
    // adds point j to tour at position i
    tour[i] = j;
    visited[j] = true;    
  }
  
  private int closestUnvisited(int i) {
    // return closest: points[closest] is the nearest unvisited point to points[i]
    Point p = points[i];
    double px = p.x();
    double py = p.y();
    int closest = -1;
    double mind = Double.MAX_VALUE;
    for (int j = 0; j < N; j++) {
      if (!visited[j]) {

        // is point j outside the outer enclosing box for prev, ie not closer than min point?
        // this is a measure of "farther" without calling distTo
        if (points[j].x() > px + mind || points[j].x() < px - mind) continue;
        if (points[j].y() > py + mind || points[j].y() < py - mind) continue;

        double d = p.distTo(points[j]);
        if (d < mind) {
          mind = d;
          closest = j;
        }          
      }
    }
     
    return closest;
       
  }

  private void twoOpt() {
    // 2-opt local improvement - exhaustive search
    int twoOptStart = (fixed) ? 1 : 0;
    for (int i = twoOptStart; i < N; i++) {      
      for (int j = twoOptStart; j < N; j++) {
        if (i == j) continue;
        if (i < j) trySwap(i, j);
        else       trySwap(j, i);
      }
    }    
  }
  
  private void trySwap(int i, int j) {
    // does a 2-opt swap on i, j if it improveds teh objective
    // assumes i < j
    
    int ip1 = ((i + 1) < N) ? i + 1 : 0;
    int jp1 = ((j + 1) < N) ? j + 1 : 0;
    double di = points[tour[i]].distTo(points[tour[ip1]]);   // length from i
    double dj = points[tour[j]].distTo(points[tour[jp1]]);
    double newdi = points[tour[i]].distTo(points[tour[j]]);  // new length from i
    double newdj = points[tour[ip1]].distTo(points[tour[jp1]]);
    double delta = - (di + dj) + newdi + newdj;
    
    if (delta < 0 || aneal(delta)) {  
      int[] temp = new int[j-i];
      for (int k = 0; k < j-i; k++) temp[k] = tour[j-k];
      for (int k = 0; k < j-i; k++) tour[i+1+k] = temp[k];
      tourLength += delta;
    }
  }
  
  private void trySwap(int i, int j, int k) {
    // 3-opt swap om i, j, k
    // assumes i < j < k
    
    int ip1 = ((i + 1) < N) ? i + 1 : 0;
    int jp1 = ((j + 1) < N) ? j + 1 : 0;
    int kp1 = ((k + 1) < N) ? k + 1 : 0;
    double di = points[tour[i]].distTo(points[tour[ip1]]);   // current segment lengths
    double dj = points[tour[j]].distTo(points[tour[jp1]]); 
    double dk = points[tour[k]].distTo(points[tour[kp1]]);          
    double newdi1 = points[tour[i]].distTo(points[tour[jp1]]);  // new segment lengths
    double newdj1 = points[tour[j]].distTo(points[tour[kp1]]);                    
    double newdk1 = points[tour[k]].distTo(points[tour[ip1]]);
    double delta1 = - (di + dj + dk) + newdi1 + newdj1 + newdk1;
    double newdi2 = points[tour[i]].distTo(points[tour[jp1]]);  // alternate new segment lengths
    double newdj2 = points[tour[k]].distTo(points[tour[j]]);                    
    double newdk2 = points[tour[kp1]].distTo(points[tour[ip1]]);
    double delta2 = - (di + dj + dk) + newdi2 + newdj2 + newdk2; 

    if (delta1 < 0 || delta2 < 0) {
      if (delta1 < delta2) {
        int[] tempi = new int[j-i];
        int[] tempj = new int[k-j];
        for (int l = 0; l < k-j; l++) tempj[l] = tour[j+1+l];
        for (int l = 0; l < j-i; l++) tempi[l] = tour[i+1+l];
        for (int l = 0; l < k-j; l++) tour[i+1+l] = tempj[l];
        for (int l = 0; l < j-i; l++) tour[i+k-j+1+l] = tempi[l];             
        tourLength += delta1;
           
      } else {
         int[] tempi = new int[j-i];
         int[] tempj = new int[k-j];
         for (int l = 0; l < k-j; l++) tempj[l] = tour[j+1+l];
         for (int l = 0; l < j-i; l++) tempi[l] = tour[j-l];
         for (int l = 0; l < k-j; l++) tour[i+1+l] = tempj[l];
         for (int l = 0; l < j-i; l++) tour[i+k-j+1+l] = tempi[l];  
         tourLength += delta2;             
      }
    }    
  }
  
  private boolean aneal(double delta) {
    if (!anneal || temperature < temperatureMin) return false;
    boolean accept = (rng.nextDouble() < Math.exp(-delta / temperature));
    if (accept && temperatureCount++ > temperatureChange) {
      temperature *= temperatureDelta;
      temperatureCount = 0;
    }
    return accept;
  }
  
  private void randomTwoOpt(int multiplier) {
    for (int k = 0; k < multiplier * N; k++) {
      
      int i = rng.nextInt(N); 
      int j = rng.nextInt(N);  

      if (fixed) if (i == 0 || j == 0) continue;
      if (i == j) continue;
      
      if (i < j) trySwap(i, j);
      else       trySwap(j, i);
    }
  }
  
  private void randomThreeOpt(int multiplier) {
    for (int a = 0; a < multiplier * N; a++) {
      
      int i = rng.nextInt(N); 
      int j = rng.nextInt(N);
      int k = rng.nextInt(N);  

      if (fixed) if (i == 0 || j == 0 || k == 0) continue;
      if (i == j || j == k || i == k) continue;
      
      if      (i < j && j < k) trySwap(i, j, k);
      else if (i < k && k < j) trySwap(i, k, j);
      else if (j < i && i < k) trySwap(j, i, k);
      else if (j < k && k < i) trySwap(j, k, i);
      else if (k < i && i < j) trySwap(k, i, j);
      else                     trySwap(k, j, i);
    }
  }

  private void threeOpt() {
    // 3-opt local improvement - exhaustive search ... SLOW!
    int twoOptStart = (fixed) ? 1 : 0;
    for (int i = twoOptStart; i < N; i++) {      
      for (int j = twoOptStart; j < N; j++) {
        for (int k = twoOptStart; k < N; k++) {
          
          if (i == j || j == k || i == k) continue;
          
          if      (i < j && j < k) trySwap(i, j, k);
          else if (i < k && k < j) trySwap(i, k, j);
          else if (j < i && i < k) trySwap(j, i, k);
          else if (j < k && k < i) trySwap(j, k, i);
          else if (k < i && i < j) trySwap(k, i, j);
          else                     trySwap(k, j, i);
        } 
      }
    }   
  }

  public boolean feasible() {
    // is the tour feasible, eg is each point in the tour exactly one
    // sum of tour points should equal sum of N-1 integers
    int count = 0;
    for (int i = 0; i < N; i++) count += tour[i];
    return (count == ((N) * (N-1) / 2));
  }
  
  public int[] tour() { return tour; }
  public void setMap(int[] m) { map = m; }
  public int[] mappedTour() {
    int[] mappedTour = new int[N];
    for (int j = 0; j < N; j++) {
      mappedTour[j] = map[tour[j]];
    } 
    return mappedTour;   
  }
  
  public double length() {
    double l = 0;
    for (int i = 0; i < tour.length-1; i++) {
      l += points[tour[i]].distTo(points[tour[i+1]]);
    }
    return l + points[tour[N-1]].distTo(points[tour[0]]);
  }
  
  public double maxLeg() {
    double max = Double.MIN_VALUE;
    for (int i = 0; i < tour.length-1; i++) {
      double l = points[tour[i]].distTo(points[tour[i+1]]);
      if (l > max) max = l;
    }
    double l = points[tour[N-1]].distTo(points[tour[0]]);
    if (l > max) max = l;
    return max;
  }
  
  public double minLeg() {
    double min = Double.MAX_VALUE;
    for (int i = 0; i < tour.length-1; i++) {
      double l = points[tour[i]].distTo(points[tour[i+1]]);
      if (l < min) min = l;
    }
    double l = points[tour[N-1]].distTo(points[tour[0]]);
    if (l < min) min = l;
    return min;
  }
  
  public double meanLeg() { return length() / (N + 1); }
  
  public String annealingParams() {
    String s = "starting temperature = " + temperatureStart + "\n";
    s += "minimum temperature = " + temperatureMin + "\n";
    s += "delta = " + temperatureDelta + "\n";
    s += "step = " + temperatureChange;
    return s;
  }
  
  public double time() { return time; }
    
  public String toString() {
    String s = length() + " 0\n";      
    for (int i = 0; i < tour.length; i++) {
      s += tour[i] + " ";
    }
    return s;
  }

  private class Point {
    private double x, y;
    int point;
    private Point (double a, double b, int i) {
      x = a;
      y = b;
      point = i;
    }
    private double x() { return x; }
    private double y() { return y; }
    private int point() { return point; }
    private double distTo (Point that) {
      double tx = x - that.x();
      double ty = y - that.y();
      return Math.sqrt(tx * tx + ty * ty);
    }
    public String toString() {
      return "(" + x + ", " + y + ")";
    }
  }

}