import java.io.*;
import java.util.List;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Random;


public class Solver {
  
  private int N, M;
  private Warehouse[] warehouses; 
  private Customer[] customers;
  private int newWarehouseMode = 1;

    
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
    M = Integer.parseInt(firstLine[1]);
    
    warehouses = new Warehouse[N];
    customers = new Customer[M]; 

    for (int j = 0; j < N; j++){
      String line = lines.get(j + 1);
      String[] parts = line.split("\\s+");
      warehouses[j] = new Warehouse(j, Double.parseDouble(parts[0]), Double.parseDouble(parts[1]));
    }  
    
    for (int i = 0; i < M; i++) {
      String line = lines.get(N + 2*i + 1);
      String[] parts = line.split("\\s+");
      double demand = Double.parseDouble(parts[0]);
      line = lines.get(N + 2*i + 2);
      parts = line.split("\\s+");
      double[] costs = new double[N];
      for (int j = 0; j < N; j++) costs[j] = Double.parseDouble(parts[j]);
      customers[i] = new Customer(i, demand, costs);
    }
    
    /*
    *     a greedy solution
    */
    
    // start with all warehouses closed
    for (int j = 0; j < N; j++) warehouses[j].setClosed();
    
    // construct a greedy feasible solution
    //java.util.Arrays.sort(customers);
    while (!feasible(false) && numOpenWarehouses() < N) {
      openNewWarehouse(newWarehouseMode);
      clearWarehouses();
      assignCustomers();
      openSwap();
    }
    refine123();
    feasible(true);
    printSolution(true);

    // refine the solution by swapping closed/open warehouses
    // todo
    closedSwap();
    //closedSwap();
    
    // print the full solution
    printSolution(true);
   
  } 
  
  private void closedSwap() {
    Bag<Integer> openWarehouses = new Bag<Integer>();
    Bag<Integer> closedWarehouses = new Bag<Integer>();    
    for (int j = 0; j < N; j++) {
      if (warehouses[j].open()) openWarehouses.add(j);
      else                      closedWarehouses.add(j);
    }
    
    //String s1 = "closed: ";
    //for (int closedIndex : closedWarehouses) s1 += closedIndex + " ";
    //System.out.println(s1);
    //String s2 = "open: ";
    //for (int openIndex : openWarehouses) s2 += openIndex + " ";
    //System.out.println(s2);
        
    for (int closedIndex : closedWarehouses) {
      for (int openIndex : openWarehouses) {
        swap(warehouses[closedIndex], warehouses[openIndex]);
      }
    }
  }
  
  private void openSwap() {
    Bag<Integer> openWarehouses = new Bag<Integer>();
    for (int j = 0; j < N; j++) {
      if (warehouses[j].open()) openWarehouses.add(j);
    }
    
    //String s2 = "open: ";
    //for (int openIndex : openWarehouses) s2 += openIndex + " ";
    //System.out.println(s2);
        
    for (int openIndex1 : openWarehouses) {
      for (int openIndex2 : openWarehouses) {
        if (openIndex1 == openIndex2) continue;
        swap(warehouses[openIndex1], warehouses[openIndex2]);
      }
    }
  }
  
  private void swap(Warehouse w1, Warehouse w2) {
    // swap customers assigned to warehouses w1, w2 if it improves the solution
    int windex1 = w1.windex();
    int windex2 = w2.windex();
    
    Bag<Integer> w2Customers = new Bag<Integer>();
    for (int i = 0; i < M; i++) {
      if (customers[i].windex() == windex2) {
        w2Customers.add(customers[i].cindex());
      }
    }
    
    double sum1 = warehouses[windex1].setup();
    double sum2 = warehouses[windex2].setup();
    for (int cindex : w2Customers) {
      sum1 += customers[cindex].cost(windex1);
      sum2 += customers[cindex].cost(windex2);
    }
    //System.out.println("w1, w2, sum1, sum2 = " + windex1 + ", " + windex2 + ", " + sum1 + ", " + sum2);
    //String s = "w2 customers: ";
    //for (int cindex : w2Customers) s += cindex + " ";
    //System.out.println(s);
        
    // will the swap improve?
    if (sum1 < sum2) {

      
      //yes, make the swap
      for (int cindex : w2Customers) {
        customers[cindex].assignWarehouse(warehouses[windex1]);
      }
      warehouses[windex1].setOpen();
      warehouses[windex2].setClosed(); 
      //printSolution(true);    
    }
    
  }
  
  private void assignCustomers() {
    // assigns customers greedily to open warehouses
    // note: may not be feasible in that all customers may not be assigned
    for (int i = 0; i < M; i++) {
      Customer c = customers[i];
      MinPQ<WarehouseCost> pq = costsPQ(c);
      for (WarehouseCost wc : pq) {
        Warehouse w = warehouses[wc.windex()];
        if (w.open() && w.fits(c)) {
          c.assignWarehouse(w);
          break;
        }
      }
      
    }        
  }
  
  private void clearWarehouses() {
    // clears all customers warehouse assignments
    for (int i = 0; i < M; i++) {
      customers[i].clearWarehouse();
    }
  }
  
  
  private MinPQ<WarehouseCost> costsPQ(Customer c) {
    MinPQ<WarehouseCost> pq = new MinPQ<WarehouseCost>();
    for (int j = 0; j < N; j++) pq.insert(new WarehouseCost(j, c.cost(j)));
    return pq;    
  }
  
  private void refine123() {
    for (int i = 0; i < M; i++) if (customers[i].assigned()) refine123(customers[i]);    
  }
  
  private boolean refine123(Customer c) {
    // fits 1-cheap?
    Warehouse wc1 = warehouses[c.minWindex()];
    if (wc1.open() && wc1.fits(c)) {
      c.assignWarehouse(wc1);
      return true;
    }
    
    // fits 2-cheap?
    double min = Double.MAX_VALUE;
    Warehouse wc2 = null;
    for (int j = 0; j < N; j++) {
      if (warehouses[j] == wc1) continue;
      if (c.cost(j) < min) {
        min = c.cost(j);
        wc2 = warehouses[j];
      }
    }
    if (wc2.open() && wc2.fits(c)) {
      c.assignWarehouse(wc2);
      return true;
    }
    
    // fits 3-cheap?
    min = Double.MAX_VALUE;
    Warehouse wc3 = null;
    for (int j = 0; j < N; j++) {
      if (warehouses[j] == wc1 || warehouses[j] == wc2) continue;
      if (c.cost(j) < min) {
        min = c.cost(j);
        wc3 = warehouses[j];
      }
    }
    if (wc3.open() && wc3.fits(c)) {
      c.assignWarehouse(wc3);
      return true;
    }   
    
    return false;
     
  }
  
  private void refine1() {
    
    // refines the current solution by searching for opportunities to change 
    // customer warehouse assignments
    //openNewWarehouse(newWarehouseMode);
    double z = obj();
    for (int i = 0; i < M; i++) {
      Customer c = customers[i]; 
      if (!c.assigned()) continue;                // customer not assigned yet     
      if (c.windex() == c.minWindex()) continue;  // can't improve assignment
      MinPQ<WarehouseCost> pq = costsPQ(c);
      for (WarehouseCost wc : pq) {
        if (wc.windex() == c.windex()) continue;
        Warehouse w = warehouses[wc.windex()];
        if (w.open() && w.fits(c) && wc.cost() < c.cost(c.windex())) c.assignWarehouse(w);
        else break;
      }
    }
 
  }
  
  private int openNewWarehouse(int mode) {
    // opens the currently closed warehouse with max efficiency (bang for buck), capacity, lowest setup, etc.
    // returns the warehouses index for the newly opened warehouse
    
    // mode = 1 => max efficiency
    // mode = 2 => max capacity
    // mode = 3 => min setup cost
    
    double max1 = Double.MIN_VALUE;
    double max2 = Double.MIN_VALUE;
    double min3 = Double.MAX_VALUE;
    int w1 = -1;
    int w2 = -1;
    int w3 = -1;
    
    for (int j = 0; j < N; j++) {
      Warehouse w = warehouses[j];
      if (!w.open()) {
        if (w.efficiency() > max1) {
          max1 = w.efficiency();
          w1 = j;
        }
        if (w.capacity() > max2) {
          max2 = w.capacity();
          w2 = j;
        }
        if (w.setup() < min3) {
           min3 = w.setup();
           w3 = j;
         }
      }
    }
    
    int w = -1;
    if      (mode == 1) w = w1;
    else if (mode == 2) w = w2;
    else                w = w3;
    if (w != -1 ) warehouses[w].setOpen();
    return w;
  }
  
  private void printSolution(boolean printFull) {
    // print the solution in the specified output format        
    double obj = 0.0;
    boolean[] warehouseUsed = new boolean[N];
    String s = "";
    for (int i = 0; i < M; i++) {
      int windex = customers[i].windex();
      if (windex != -1) {
        warehouseUsed[windex] = true;
        obj += customers[i].cost(windex);;
      }
      s += windex + " ";
    }
    for (int j = 0; j < N; j++) {
      if (warehouseUsed[j]) obj += warehouses[j].setup();
    }
    System.out.println(obj + " 0");
    if (printFull) System.out.println(s);
  }
  
  private double obj() {
    // return the objective function for the current solution        
    double objc = 0.0;
    boolean[] warehouseUsed = new boolean[N];
    for (int i = 0; i < M; i++) {
      if (!customers[i].assigned()) continue;
      int windex = customers[i].windex();
      warehouseUsed[windex] = true;
      objc += customers[i].cost(windex);
    }
    double objw = 0;
    for (int j = 0; j < N; j++) {
      if (warehouseUsed[j]) objw += warehouses[j].setup();
    }
    //System.out.println("objc, objw = " + objc + ", " + objw);
    return objc + objw;
  }
  
  private boolean feasible(boolean verbose) {
    // is the current solution feasible?
    // all customers assigned to a warehouse and all warehouses under capacity
    for (int i = 0; i < M; i++) {
      if (!customers[i].assigned()) {
        if (verbose) System.out.println("customer not assigned = " + customers[i]);
        return false;
      }
    }
    for (int j = 0; j < N; j++) {
      if (warehouses[j].spare() < 0) {
        if (verbose) System.out.println("warehouse over capacity = " + j);
        return false;        
      }
    }
    return true;
  } 
  
  private int numOpenWarehouses() {
    int count = 0;
    for (int j = 0; j < N; j++) if (warehouses[j].open()) count++;
    return count;
  }

  /*
  *  the warehouse model
  */
  private class Warehouse {
    private double capacity, used, setup;
    int windex;     // warehouses[] index of this warehouse
    boolean open;
    public Warehouse (int j, double c, double w) {
      windex = j;
      capacity = c;
      setup = w;
      used = 0.0;
      open = false;
    }
    private double capacity() { return capacity; }
    private void incrementUsed(double d) { used += d; }
    private void decrementUsed(double d) { used -= d; }
    private double spare() { return capacity - used; }
    private double setup() { return setup; }
    private int windex() { return windex; }
    private double efficiency() { return setup / capacity; }
    private boolean open() { return open; }
    private void setClosed() { open = false; }
    private void setOpen() { open = true; }
    private boolean available() { return open && (capacity > used); }
    private boolean fits(Customer c) { return c.demand() <= spare(); }
    public String toString() {
      return "warehouse, capacity, setup = " + windex + ", " + capacity + ", " + setup;
    }
  }
  
  /*
  * the customer model
  */
  private class Customer implements Comparable<Customer> {
    private double demand;
    private double[] costs;
    private int cindex;       // customers[] index of this customer
    private int windex;       // warehouses[] index of assigned warehouse, -1 if unassigned
    private int minWindex;
    private boolean assigned;
    private Customer(int i, double d, double[] c) {
      windex = -1;
      assigned = false;
      cindex = i;
      demand = d;
      costs = new double[N];
      double min = Double.MAX_VALUE;
      minWindex = -1;
      for (int j = 0; j < N; j++) {
        costs[j] = c[j];
        if (costs[j] < min) {
          min = costs[j];
          minWindex = j;
        }
      }
    }
    public int compareTo(Customer that) {
      if (this.demand > that.demand) return -1;
      if (this.demand < that.demand) return 1;
      return 0;
    }
    private int cindex() { return cindex; }
    private double demand() { return demand; }
    private double cost(int j) { return costs[j]; }
    private double[] costs() { return costs; }
    private int minWindex() { return minWindex; }
    private double minCost() { return costs[minWindex]; }
    private int windex() { return windex; }
    private boolean assigned() { return assigned; }
    private void assignWarehouse(Warehouse w) { 
      if (windex != -1) clearWarehouse();
      windex = w.windex(); 
      w.incrementUsed(demand); 
      assigned = true; 
    }
    private void clearWarehouse() { 
      if (windex != -1) warehouses[windex].decrementUsed(demand);
      windex = -1;
      assigned = false;
    }
    private boolean fits(Warehouse w) {
      // can customer "fit" ie be supplied by the warehouse
      return (demand <= w.spare());      
    }
    public String toString() {
      String s = "";
      if (assigned) s = "customer, demand, warehouse, cost = " + cindex + ", " + demand + ", " + windex + ", " + costs[windex];
      else          s = "customer, demand, warehouse = " + cindex + ", " + demand + ", " + windex; 
      return s;
    }

    
  }
  
  private class WarehouseCost implements Comparable<WarehouseCost>{
    private int windex;
    private double cost;
    private WarehouseCost(int w, double c) {
      windex = w;
      cost = c;
    }
    private int windex() { return windex; }
    private double cost() { return cost; }
    public int compareTo(WarehouseCost that) {
      if (this.cost > that.cost) return 1;
      if (this.cost < that.cost) return -1;
      return 0;
    }
  }
  

  
}