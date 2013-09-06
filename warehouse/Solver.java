import java.io.*;
import java.util.List;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Random;

/**
 * The class <code>Solver</code> is an implementation of a greedy algorithm to solve the knapsack problem.
 *
*/

public class Solver {
  
  // data defining the customer model
  int M;
  double[][] cost;
  double[] demand;
  int[] warehouse;
  boolean[] open;
  
  // data defining the warehouse model
  int N;
  double[] setup;
  double[] capacity;
  double[] spare;
  
  // parameters for simulated annealing
  boolean anneal = false;
  double temperatureStart = 1000;
  double temperature = temperatureStart;
  double temperatureMin = 0.0001;
  double temperatureDelta = 1 - temperatureMin;
  double temperatureCount = 0;
  double temperatureChange = 100;
  
  // other data
  Random rng = new Random();

    
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
    N = Integer.parseInt(firstLine[0]);   // number of warehoused
    M = Integer.parseInt(firstLine[1]);   // number of customers
    
    // data defining the warehouse model
    setup = new double[N];
    capacity = new double[N];
    spare = new double[N];
    open = new boolean[N];

    for (int j = 0; j < N; j++){
      String line = lines.get(j + 1);
      String[] parts = line.split("\\s+");
      capacity[j] = Double.parseDouble(parts[0]);
      setup[j] = Double.parseDouble(parts[1]);
      spare[j] = capacity[j];
    } 
    
    // data defining the customer model
    cost = new double[M][N];
    demand = new double[M];
    warehouse = new int[M]; 
    
    for (int i = 0; i < M; i++) {
      String line = lines.get(N + 2*i + 1);
      String[] parts = line.split("\\s+");
      demand[i] = Double.parseDouble(parts[0]);
      line = lines.get(N + 2*i + 2);
      parts = line.split("\\s+");
      for (int j = 0; j < N; j++) cost[i][j] = Double.parseDouble(parts[j]);
      warehouse[i] = -1;
    }
    

    closeAll();

    greedy(7);
    refine2();
    feasible();
    printSolution(false);   

    refine2();
    
    feasible(); 
    printSolution();   
  } 
  
  private void greedy() { 
    greedy(openWarehouses()); 
  }
  
  private void greedy(int n) { 
    openNextMaxWarehouses(n); 
    greedy();
  }
  
  private void greedy(Bag<Integer> initialOpenWarehouses) {
    
    // greedy: assign each customer to the open warehouse that can most cheaply supply it
    // handle nuances about warehouse capacity, etc.
    // start with some number of warehouses open
    
    for (int j : initialOpenWarehouses) open(j);
    assignToBestOpen();
    
    for (int i = 0; i < M; i++) {
      while(warehouse[i] == -1) {
        // a customer is not served, eg did not fit the open warehouses
        // open a new warehouse and reassign
        if (openNextMaxWarehouses(1) == 0) break;
        assignToBestOpen();
      }
    }

    // close any warehouses that were opened but not used
    for (int j = 0; j < N; j++) if (spare[j] == capacity[j]) close(j);
  }
  
  private void refine1() {
    // refine the solution by trying the nest 'best' warehouses
    Bag<Integer> openWarehouses = openWarehouses();
    Bag<Integer> bestWarehouses = new Bag<Integer>();
    openNextMinWarehouses(1);  
    double objective = obj();
    for (int j1 : openWarehouses) {
      close(j1);
      assignToBestOpen();
      double z = obj();
      if (z < objective) {
        objective = z;
        bestWarehouses = openWarehouses();
      }
      open(j1);
    }
    
    closeAll();
    for (int j : bestWarehouses) open(j);
    assignToBestOpen();
  }
  
  private void refine2() {
    // refine the solution by trying the nest 'best' warehouses
    Bag<Integer> openWarehouses = openWarehouses();
    debug(openString());
    
    Bag<Integer> bestWarehouses = new Bag<Integer>();
    openNextMinWarehouses(2);  
    double objective = obj();
    for (int j1 : openWarehouses) {
      close(j1);
      for (int j2 : openWarehouses) {
        if (j2 == j1) continue;
        close(j2);
        assignToBestOpen();
        double z = obj();
        if (z < objective) {
          objective = z;
          bestWarehouses = openWarehouses();
        }
        open(j2);        
      }
      open(j1);
    }
    
    closeAll();
    for (int j : bestWarehouses) open(j);
    debug(openString());
    assignToBestOpen();
  }
  
  
  /*
  * customer and warehouse swaps
  */
  
  private void doCustomerSwaps() {
    // try to improve the objective by swapping pairs of customers among the open warehouses
    for (int i1 = 0; i1 < N; i1++) {
      for (int i2 = 0; i2 < N; i2++) {
        tryCustomerSwap(i1, i2);
      }
    }
  }
  
  private void doRandomCustomerSwaps(int multiplier) {
    // try to improve score with random customer pair swaps
    for (int k = 0; k < multiplier * N; k++) tryCustomerSwap(rng.nextInt(N), rng.nextInt(N));
  }
  
  private void tryCustomerSwap(int i1, int i2) {
    // swaps warehouse assignment for customers i1 and i2 if the objective improves
    int j1 = warehouse[i1];
    int j2 = warehouse[i2];
    if (j1 == j2) return;   // i1 and i2 must be in different warehouses (and thus i1 != i2)
    
    double s1 = spare[j1] + demand[i1] - demand[i2];
    double s2 = spare[j2] + demand[i2] - demand[i1];
    if (s1 < 0 || s2 < 0) return;   // one of the swapees won't fit
    
    double d =  - cost[i1][j1] - cost[i2][j2] + cost[i1][j2] + cost[i2][j1];

    if (d < 0 || anneal(d)) {
      // make the swap
      assign(i1, j2);
      assign(i2, j1);
    }
    
  }
  
  private void doWarehouseSwaps() {
    // tries to improve objective by swapping the customers served by an open warehouse with another warehouse
    for (int j1 = 0; j1 < N; j1++) {
      if (!open[j1]) continue;
      for (int j2 = 0; j2 < N; j2++) {
        tryWarehouseSwap(j1, j2);
      }
    }
  }
  
  private void doRandomWarehouseSwaps(int multiplier) {
    // try to improve score with random customer pair swaps
    for (int k = 0; k < multiplier * N; k++) tryWarehouseSwap(rng.nextInt(N), rng.nextInt(N));
  }
  
  private void tryWarehouseSwap(int j1, int j2) {
    // swaps customers among warehouses
    // assumes j1 & j2 are open
    if (j1 == j2) return;
    
    double cost1 = 0;
    double demand1 = 0;
    Bag<Integer> bag1 = new Bag<Integer>();
    for (int i = 0; i < M; i++) {
      if (warehouse[i] == j1) {
        cost1 += cost[i][j1];
        demand1 += demand[i];
        bag1.add(i);
      }
    }
    
    double cost2 = 0;
    double demand2 = 0;
    Bag<Integer> bag2 = new Bag<Integer>();
    for (int i = 0; i < M; i++) {
      if (warehouse[i] == j2) {
        cost2 += cost[i][j2];
        demand2 += demand[i];
        bag2.add(i);
      }
    }

    // will the swap fit both warehouses?
    if (demand1 > capacity[j2] || demand2 > capacity[j1]) return;
    
    double cost12 = 0;
    for (int i : bag1) cost12 += cost[i][j2];
    double cost21 = 0;
    for (int i : bag2) cost21 += cost[i][j1];
    
    double delta = cost12 + cost21 - cost1 - cost2;
    if (delta < 0 || anneal(delta)) {
      
      // make the swap
      for (int i : bag1) {
        clear(i);
        assign(i, j2);
      }
      for (int i : bag2) {
        clear(i);
        assign(i, j1);
      }
      
      if(!open[j2]) {
        open(j2);
        close(j1);
      }
    }
                  
  }
  
  private void closedSwap() {
    Bag<Integer> openWarehouses = new Bag<Integer>();
    Bag<Integer> closedWarehouses = new Bag<Integer>();    
    for (int j = 0; j < N; j++) {
      if (open[j]) openWarehouses.add(j);
      else         closedWarehouses.add(j);
    }
    
    for (int j1 : closedWarehouses) {
      for (int j2 : openWarehouses) {
        swapWarehouses(j1, j2);
      }
    }
  }
  
  private void openSwap() {
    Bag<Integer> openWarehouses = new Bag<Integer>();
    for (int j = 0; j < N; j++) {
      if (open[j]) openWarehouses.add(j);
    }
        
    for (int j1 : openWarehouses) {
      for (int j2 : openWarehouses) {
        if (j1 == j2) continue;
        swapWarehouses(j1, j2);
      }
    }
  }
  
  private void swapWarehouses(int j1, int j2) {
    // swap customers assigned to warehouses j1, j2 if it improves the solution
    // assumes j2 is open
    boolean openj1 = open[j1];
    boolean openj2 = open[j2];
    
    Bag<Integer> j1Customers = assignedCustomers(j1);
    Bag<Integer> j2Customers = assignedCustomers(j2);
    
    double sum11 = 0;
    double sum12 = 0;
    if (openj1) {
      sum11 += setup[j1];
      sum12 += setup[j2];
      for (int i : j1Customers) {
        sum11 += cost[i][j1];
        sum12 += cost[i][j2];
      } 
    } 
        
    double sum21 = 0;
    double sum22 = 0;
    if (openj2) {
      sum21 += setup[j1];
      sum22 += setup[j2];
      for (int i : j2Customers) {
        sum21 += cost[i][j1];
        sum22 += cost[i][j2];
      }      
    }
        
    // will the swap improve?
    if (sum21 + sum12 < sum22 + sum11) {
      
      //yes, make the swap
      if (openj2) {
        // j2 customers to j1
        open(j1);
        for (int i : j2Customers) assign(i, j1);
      } else {
        close(j1);
      }

      if (openj1) {
        // j1 customers to j2
        open(j2);
        for (int i : j1Customers) assign(i, j2);        
      } else {
        close(j2); 
      }
    }
    
  }
  
  /*
  * assigning and fitting
  */

  private boolean tryFit(int i) {
    // tries to fit customer i into the open warehouses by shuffling service assignments
    boolean fit = false;
    
    
      
    return fit;
  }
  
  private void assignToBestOpen() {
    // clears all warehouse assignments and reassigns customers to best (cheapest) open warehouse
    // eg, used after a new warehouse is open to try to improve overall assignments
    for (int i = M-1; i >= 0; i--) clear(i);
    for (int i = M-1; i >= 0; i--) assignToBestOpen(i);
  }
  
  private void assignToBestOpen(int i) {
    // assigns customer i to the best open warehouse
    int j = bestOpenFit(i);
    if (j != -1) assign(i, j);
  }
  
  /* opening warehouses
  */
  
  private int openNextMaxWarehouses(int n) {
    // opens the n closed warehouses with max bang for buck
    int count = 0;
    for (int k = 0; k < n; k++) if (openNextMaxWarehouse() != -1) count++;
    return count;
  }
  
  private int openNextMaxWarehouse() {
    // opens one closed warehouse with max bang for buck
    double max = Double.MIN_VALUE;
    int maxj = -1;
    for (int j = 0; j < N; j++) {
      if (open[j]) continue;
      double e = bang(j);
      if (e > max) {
        max = e;
        maxj = j;
      }
    }
    if (maxj != -1) open(maxj);
    return maxj;
  }
  
  private int openNextMinWarehouses(int n) {
    // opens the n closed warehouses with min cost (expected, total serve, etc.)
    int count = 0;
    for (int k = 0; k < n; k++) if (openNextMinWarehouse() != -1) count++;
    return count;
  }
  
  private int openNextMinWarehouse() {
    // opens the n closed warehouse with min cost (expected, total serve, etc.)
    double min = Double.MAX_VALUE;
    int minj = -1;
    for (int j = 0; j < N; j++) {
      if (open[j]) continue;
      double e = setup[j];
      if (e < min) {
        min = e;
        minj = j;
      }
    }
    if (minj != -1) open(minj);
    return minj;
  }
  
  private int openNextRandomWarehouses(int n) {
    // opens the n closed warehouses with min cost (expected, total serve, etc.)
    int count = 0;
    for (int k = 0; k < n; k++) if (openNextRandomWarehouse() != -1) count++;
    return count;    
  }
  
  private int openNextRandomWarehouse() {
    int numOpen = 0;
    for (int j = 0; j < N; j++) if (open[j]) numOpen++;
    if (numOpen == N) return -1;
    int k = rng.nextInt(N - numOpen);
    // open the kth unopen warehouse encountered
    int l = 0;
    for (int j = 0; j < N; j++) {
      if (!open[j]) {
        if (l++ == k) {
          open(j);
          return j;
        }
      }
    }
    
    return -1;
    
  }
  
  private void openAllWarehouses() {
    // opens all warehouses (currently open warehouses are not touched)
    for (int j = 0; j < N; j++) if (!open[j]) open(j);
  }
  
  private double bang(int j) { return capacity[j] / setup[j]; }
  
  private double expected(int j) {
    // returns the expected cost to open and operate warehouse j
    double dbar = 0;
    double tbar = 0;
    for (int i = 0; i < M; i++) {
      dbar += demand[i];
      tbar += cost[i][j];
    }
    dbar /= M;
    tbar /= M;
    return setup[j] + capacity[j] / dbar + tbar;
  }
  
  private Bag<Integer> openWarehouses() {
    // returns bag containing all open warehouses
    Bag<Integer> b = new Bag<Integer>();
    for (int j = 0; j < N; j++) if (open[j]) b.add(j);
    return b;
  }
  
  private Bag<Integer> closedWarehouses() {
    // returns bag containing all closed warehouses
    Bag<Integer> b = new Bag<Integer>();
    for (int j = 0; j < N; j++) if (!open[j]) b.add(j);
    return b;
  }
  
  private Bag<Integer> assignedCustomers(int j) {
    // returns a bag containing all the customers assigned to warehouse j
    Bag<Integer> b = new Bag<Integer>();
    for (int i = 0; i < M; i++) if (warehouse[i] == j) b.add(i);
    return b;    
  }
  
  private double totalDemand() {
    double sum = 0;
    for (int i = 0; i < M; i++) sum += demand[i];
    return sum;
  }
  
  private double totalOpenCapacity() {
    double sum = 0;
    for (int j: openWarehouses()) sum += capacity[j];
    return sum;
  }
  
  /*
  * simulated annealing
  */
  
  private boolean anneal(double delta) {
    if (!anneal || temperature < temperatureMin) return false;
    boolean accept = (rng.nextDouble() < Math.exp(-delta / temperature));
    if (accept && temperatureCount++ > temperatureChange) {
      temperature *= temperatureDelta;
      temperatureCount = 0;
    }
    return accept;
  }
  
  /*
  * helpers and misc
  */
  
  private double serveAllCost(int j) {
    // returns the total cost if warehouse j were to serve all customers (ignoring capacity constraints)
    double sum = 0;
    for (int i = 0; i < M; i++) sum += cost[i][j];
    return setup[j] + sum;
  }
  
  private boolean feasible() {
    // is the current solution feasible?
    // all customers assigned to a warehouse and all warehouses under capacity
    for (int i = 0; i < M; i++) {
      if (warehouse[i] == -1) {
        System.out.println("customer not assigned = " + i);
        return false;
      }
    }
    for (int j = 0; j < N; j++) {
      if (spare[j] < 0) {
        System.out.println("warehouse over capacity = " + j);
        return false;        
      }
    }
    return true;
  }

  private boolean fits(int i, int j) {
    // can customer i be serviced by warehouse j?
    if (!open[j]) return false;
    return demand[i] <= spare[j]; 
  }
  
  private void assign(int i, int j) {
    // assign customer i to be serviced by warehouse j
    clear(i);
    warehouse[i] = j;
    spare[j] -= demand[i];
  }
  
  private void clear(int i) {
    // clear customer i warehouse assignment
    if (warehouse[i] != -1) {
      spare[warehouse[i]] += demand[i];
      warehouse[i] = -1;
    }
  }
  
  private void open(int j) {
    // sets warehouse j to open with full capacity (and no assigned customers)
    close(j);
    open[j] = true;
    spare[j] = capacity[j];
  }
  
  private void close(int j) {
    // sets warehouse j to closed; clears assigned customers
    for (int i = 0; i < M; i++) if (warehouse[i] == j) clear(i);
    open[j] = false;
    spare[j] = capacity[j];
  }
  
  private void closeAll() {
    // close all the warehouses and clear all the customers
    for (int j = 0; j < N; j++) close(j);
    for (int i = 0; i < M; i++) clear(i);
  }
  
  private int bestOpenFit(int i) {
    // returns the open warehouse with spare capacity that can fit customer i
    // returns -1 if not open warehouse fits customer
    double min = Double.MAX_VALUE;
    int minj = -1;
    for (int j = 0; j < N; j++) {
      if (fits(i, j) && cost[i][j] < min) {
        min = cost[i][j];
        minj = j;
      }
    }
    return minj;
  }
  
  private void printSolution(boolean verbose) {
    // print the solution in the specified output format        
    double obj = 0.0;
    boolean[] warehouseUsed = new boolean[N];
    String s = "";
    for (int i = 0; i < M; i++) {
      int j = warehouse[i];
      if (j != -1) {
        warehouseUsed[j] = true;
        obj += cost[i][j];
        s += j + " ";        
      } else {
        s += "-1 ";
      }

    }
    for (int j = 0; j < N; j++) {
      if (warehouseUsed[j]) obj += setup[j];
    }
    System.out.println(obj + " 0");
    if (verbose) System.out.println(s);
  }
  
  private void printSolution() { printSolution(true); }
  
  private void debug(String s) { System.out.println(s); }
  
  private String openString() {
    String s = "";
    for (int j = 0; j < N; j++) if (open[j]) s += j + " ";
    return s;
  }
  
  private double obj() {
    // return the objective function for the current solution        
    double objc = 0.0;
    boolean[] warehouseUsed = new boolean[N];
    for (int i = 0; i < M; i++) {
      int j = warehouse[i];
      warehouseUsed[j] = true;
      objc += cost[i][j];
    }
    double objw = 0;
    for (int j = 0; j < N; j++) {
      if (warehouseUsed[j]) objw += setup[j];
    }
    //System.out.println("objc, objw = " + objc + ", " + objw);
    return objc + objw;
  }
  
  /*
  * private classes
  */
  private class WarehouseCost implements Comparable<WarehouseCost>{
    private int warehouse;
    private double cost;
    private WarehouseCost(int w, double c) {
      warehouse = w;
      cost = c;
    }
    private int warehouse() { return warehouse; }
    private double cost() { return cost; }
    public int compareTo(WarehouseCost that) {
      if (this.cost > that.cost) return 1;
      if (this.cost < that.cost) return -1;
      return 0;
    }
  }
  

  
}