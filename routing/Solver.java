import java.io.*;
import java.util.List;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Random;

public class Solver {
  
  int N, V, K;
  Customer[] customers;
  Vehicle[] vehicles; 
  double totalRouteLength;      // the objective function
  boolean[] visited;            // has a point been visited yet on the tour? 
  int numVisited;
  TspSolver[] tsps;
  double capacity;
  Random rng;
  
  // parameters for simulated annealing
  boolean tspAnneal = false;
  boolean anneal = false;
  double temperatureStart = 10;
  double temperature = temperatureStart;
  double temperatureMin = 0.01;
  double temperatureDelta = 1 - temperatureMin;
  double temperatureCount = 0;
  double temperatureChange = 1000;
    
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
    
  public Solver(String[] args) throws IOException {
    String fileName = null;
        
    // get the temp file name
    for(String arg : args){
      if(arg.startsWith("-file=")) {
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
      while (( line = input.readLine()) != null) {
        lines.add(line);
      }
    }
    finally {
      input.close();
    }
              
    // parse the data in the file
    String[] firstLine = lines.get(0).split("\\s+");
    N = Integer.parseInt(firstLine[0]);
    V = Integer.parseInt(firstLine[1]);
    capacity = Double.parseDouble(firstLine[2]);
    customers = new Customer[N];
    vehicles = new Vehicle[V];
    for (int i = 0; i < V; i++) vehicles[i] = new Vehicle(i, capacity);

    for(int i = 0; i < N; i++) {
      String line = lines.get(i+1);
      String[] parts = line.split("\\s+");
      double demand = Double.parseDouble(parts[0]);
      Point p = new Point(Double.parseDouble(parts[1]), Double.parseDouble(parts[2]));
      customers[i] = new Customer(i, demand, p);
    }
    
    rng = new Random();
    
     
    
    /*
    * solve it
    */
/*    
    double min = Double.MAX_VALUE;
    double minx = -1;
    for (double x = 0; x < 3.1415927; x += 0.1) {
      for (int i = 0; i < V; i++) vehicles[i].setInactive();
      solverB(x);
      double trl = totalRoutesLength();
      if (trl < min) {
        min = trl;
        minx = x; 
        System.out.println("minx, trl = " + minx + ", " + trl);
      }
    }
    
    System.out.println("start big loop");
*/    
    Stopwatch timer = new Stopwatch();
    double minx = 0.1;
    double bestRoutesLength = Double.MAX_VALUE;
    for (int a = 0; a < 10; a++) {
      for (int i = 0; i < V; i++) vehicles[i].setInactive();
      solverB(minx);
      double trl = totalRoutesLength();
      System.out.println("\n" + a + " " + trl + " " + timer.elapsedTime());
      printForPlot();
      
      doRandomCustomerSwaps(100000);
      printForPlot(); 
      System.out.println(timer.elapsedTime());
      
      doRandomChainReplacements(100000);  
      printForPlot(); 
      System.out.println(timer.elapsedTime());
          
      doAdjacentRouteSwaps();
      printForPlot();
      System.out.println(timer.elapsedTime());
      
      if (!feasibleRoute(false)) continue;
      double z = Double.MAX_VALUE;
      for (int i = 0; i < 10; i++) {
        doRandomCustomerSwaps(100); 
        doRandomChainReplacements(100);
        doAdjacentRouteSwaps();
        doRouteMoves();
        printForPlot(); 
        System.out.println(timer.elapsedTime());
          
        double ztrl = totalRoutesLength();
        if (ztrl == z) break;   // no improvement, exit loop
        z = ztrl;
      }
      
      trl = totalRoutesLength();
      System.out.println(trl);
      if (trl == bestRoutesLength) break;
      if ( trl < bestRoutesLength) {
        // found a better solution, save it
        bestRoutesLength = trl;
        printSolution(true);
        
      }
    }
      
    feasibleRoute(true);
    printSolution(true);

  }
  
  private void solverA() { 
    // based on the nearest neighbor greedy, opening and filling one 
    
    // the greedy solution
    Bag<Integer> customerBag = new Bag<Integer>();
    for (int j = 1; j < N; j++) customerBag.add(j);
    visited = new boolean[N];   // has the customer been visited on some route?
    numVisited = 1;
    for (int i = 0; i < V; i++) {
      vehicles[i].setActive();
      int[] route = new int[N+1];
      nnGreedyRoute(i, route, customerBag);
      vehicles[i].setRoute(route);
      if (numVisited == N) break;
    } 
    
    // fix unvisited
    while (numVisited < N) {
      // some customer(s) not assigned to a route
      if (!addUnassignedToRoute()) break;
    } 
    
  }
  
  private void solverB(double thetaOffset) {
    // based on various clustering ideas (k-means, angle sweep)
    
    // compute the clusters
    //K = calcK();
    //K = V;
    //int[] clusters = kMeansClusters(K);
    int[] clusters = angleSweepClusters(thetaOffset);
    
    
    // assign customers in cluster k to vehicle k
    visited = new boolean[N];   // has the customer been visited on some route?
    numVisited = 1;
    for (int k = 0; k < K; k++) {
      Vehicle v = vehicles[k];
      v.setActive();
      Bag<Integer> customerBag = new Bag<Integer>();
      for (int j = 1; j < N; j++) if (clusters[j] == k) customerBag.add(j);
      int[] route = new int[N+1];
      nnGreedyRoute(v.vindex(), route, customerBag);
    }
    
    
    doTsp();
    feasibleRoute(true);
    printForPlot();
    
    Bag<Integer> unBag = new Bag<Integer>();
    for (int j = 0; j < N; j++) if (!customers[j].assigned()) unBag.add(j);
    addUnassignedToRoute(unBag);
    doTsp();
    feasibleRoute(true);
    printForPlot();
    
    for (int j = 0; j < N; j++) if (!customers[j].assigned()) addToNearestRadialRoute(j);
    doTsp();
    feasibleRoute(true);
    printForPlot();
          
    // fix any remaining unvisited - just cram them in somewhere
    while (numVisited < N) {
      // some customer(s) not assigned to a route, force them in
      if (!addUnassignedToRoute()) break;
    } 
    
    doTsp();
    feasibleRoute(true);
    printForPlot();
    
    doRouteMoves();
    feasibleRoute(true);
    printForPlot();
    
    doAdjacentRouteSwaps();
    feasibleRoute(true);
    printForPlot();
    
    doRandomCustomerSwaps(1000);
    feasibleRoute(true);
    printForPlot();
    
    //doRouteMoves();
    //printSolution(false);

  }
  
  private void solverC(double thetaOffset) {
    
    int[] clusters = angleSweepClusters(thetaOffset);
    
    // assign customers in cluster k to vehicle k
    visited = new boolean[N];   // has the customer been visited on some route?
    numVisited = 1;
    for (int k = 0; k < K; k++) {
      Vehicle v = vehicles[k];
      v.setActive();
      Bag<Integer> customerBag = new Bag<Integer>();
      customerBag.add(0);
      for (int j = 1; j < N; j++) if (clusters[j] == k) {
        v.add(j);
        customerBag.add(j);
      }
      v.setNumCustomers(customerBag.size());
      v.setRoute(doClusterTsp(customerBag));
      numVisited += v.numCustomers() - 1;
    }
    
    // fix unvisited
    while (numVisited < N) {
      // some customer(s) not assigned to a route
      if (!addUnassignedToRoute()) break;
    } 
    
    printForPlot();
    doTsp();
    printForPlot();
    doRouteMoves();
    printForPlot();
    doAdjacentRouteSwaps();
    printForPlot();
    //doRouteMoves();
    //printSolution(false);
    
  }
  
  private boolean anneal(double delta) {
    if (!anneal || temperature < temperatureMin) return false;
    boolean accept = (rng.nextDouble() < Math.exp(-delta / temperature));
    if (accept && temperatureCount++ > temperatureChange) {
      temperature *= temperatureDelta;
      temperatureCount = 0;
    }
    return accept;
  }
  
  private int[] doClusterTsp(Bag<Integer> bag) {
    // runs a tsp on the customers contained in the bag and returns the tsp route
    // 
    double[] x = new double[bag.size()];
    double[] y = new double[bag.size()];
    int[] map = new int[bag.size()];
    int j = 0;
    for (int i : bag) {
      x[j] = customers[i].location().x();
      y[j] = customers[i].location().x();
      map[j] = i;
    }
    TspSolver tsp = new TspSolver(x, y, false, true, true, tspAnneal);
    return tsp.mappedTour(map);
  }
  
  private void doTsp() {
    // send each vehicle's customers to the tsp solver
    // and update the route if improved

    double sum = 0;
    tsps = new TspSolver[V];
    for (int i = 0; i < V; i++) {
      if (!vehicles[i].active()) continue;
      Vehicle v = vehicles[i];
      int Nv = v.numCustomers();
      double[] x = new double[Nv];
      double[] y = new double[Nv]; 
      int[] vindexMap = new int[Nv];     
      for (int j = 0; j < Nv; j++) {
        vindexMap[j] = v.route(j);
        x[j] = customers[v.route(j)].location().x();
        y[j] = customers[v.route(j)].location().y();
      }
      tsps[i] = new TspSolver(x, y, false, true, true, tspAnneal);
      if (tsps[i].length() < vehicles[i].length()) {
        int[] mappedRoute = new int[Nv];
        int[] tspRoute = tsps[i].tour();
        for (int j = 0; j < Nv; j++) {
          mappedRoute[j] = vindexMap[tspRoute[j]];
        }
        v.setRoute(mappedRoute);
      }
    }
    
  }
  
  private void doRouteMoves() {
    // tries to move customers to other (active vehicle) routes to improve the total solution
    // note: other route needs to have spare capacity for the move
    
    for (int j = 1; j < N; j++) {
      Customer cj = customers[j];
      if (!cj.assigned()) continue;   // skip unassigned customers
      double dj = cj.demand();
      int vj = cj.vindex();
      
      // set up vehicle n (customer's current vehicle) by removing customer
      Vehicle vn = vehicles[vj];
      int Nvn = vn.numCustomers();
      double[] xn = new double[Nvn-1];
      double[] yn = new double[Nvn-1]; 
      int[] vnMap = new int[Nvn-1];
      int k = 0;
      for (int i = 0; i < Nvn; i++) {
        if (vn.route(i) == j) continue;
        Customer c = customers[vn.route(i)];
        xn[k] = c.location().x();
        yn[k] = c.location().y();
        vnMap[k] = vn.route(i);
        k++;
      }
      TspSolver tspn = new TspSolver(xn, yn, false, true, true, tspAnneal);
      
      boolean swapped = false;
      
      for (int i = 0; i < V; i++) {
        if (i == vj || !vehicles[i].active()) continue;
        if (vehicles[i].spare() <= dj) continue;
        
        // set up vehicle m (customer's possible new vehicle)
        Vehicle vm = vehicles[i];
        int Nvm = vm.numCustomers();
        double[] xm = new double[Nvm+1];
        double[] ym = new double[Nvm+1]; 
        int[] vmMap = new int[Nvm+1];
        for (int jj = 0; jj < Nvm; jj++) {
          int kk = vm.route(jj);
          xm[jj] = customers[kk].location().x();
          ym[jj] = customers[kk].location().y();
          vmMap[jj] = kk;
        }
        xm[Nvm] = cj.location().x();
        ym[Nvm] = cj.location().y();
        vmMap[Nvm] = j;
        TspSolver tspm = new TspSolver(xm, ym, false, true, true, tspAnneal);
 
        if (tspn.length() + tspm.length() < vn.length() + vm.length()) {
          // improved the solution, make the swap

          swapped = true; 
          
          vn.decrement(cj.demand());
          vm.increment(cj.demand());
          cj.clearVehicle();
          cj.setVindex(vm.vindex());
          int[] mappedRoute = new int[Nvn-1];
          int[] tspRoute = tspn.tour();
          for (int jj = 0; jj < Nvn-1; jj++) mappedRoute[jj] = vnMap[tspRoute[jj]];
          vn.setRoute(mappedRoute);
          vn.setNumCustomers(Nvn - 1);

          vm.setRoute(tspm.mappedTour(vmMap));
          vm.setNumCustomers(Nvm + 1);
          
          
          break;
        }
        
        if (swapped) break;
        
      }
      
    }
  }
  
  private void doAdjacentRouteSwaps() {
    // searches adjacent routes for improving swaps
    
    for (int k = 0; k < K; k++) {
      // try swaps between k and k+1
      int kp1 = (k + 1 < K) ? k + 1 : 0;
      int km1 = (k - 1 >= 0) ? k - 1 : K - 1;
      doRoutesSwap(k, kp1);
      doRoutesSwap(km1, k);
      doRoutesSwap(km1, kp1);
    }
    
  }
  
  private boolean doRoutesSwap(int v1, int v2) {
    // tries to swap customers among vehicles v1 and v2
    // get the customers for vehicles k and k+1
    Bag<Integer> v1Bag = new Bag<Integer>();
    Bag<Integer> v2Bag = new Bag<Integer>();

    for (int j = 1; j < N; j++) {
      if (customers[j].vindex() == v1) v1Bag.add(j);
      if (customers[j].vindex() == v2) v2Bag.add(j);
    }

    boolean swapped = false;
    for (int c1 : v1Bag) {
      for (int c2 : v2Bag) {
        swapped = trySwap(v1, c1, v2, c2);
        if (swapped) break;
      }
      if (swapped) break;
    }
    return swapped;
  }
  
  private void doExhaustiveRouteSwaps() {
    // swap some customers among routes to try to improv the score
    // exhaustively search pairs of routes
    
    for (int k = 0; k < K; k++) {
      for (int l = 0; l < K; l++) {
        if (l == k) continue;
        doRoutesSwap(k, l);
      }
    }
  }
  
  private void doRandomCustomerSwaps(int multiplier) {
    // randomly select customers to try to swap
    
    for (int k = 0; k < multiplier * N; k++) {
      int n = rng.nextInt(N); 
      int m = rng.nextInt(N);  
      if (n == 0 || m == 0 || n == m) continue;
      if (customers[n].vindex() == customers[m].vindex()) continue;
      trySwap(customers[n].vindex(), n, customers[m].vindex(), m);
    }
  }
  
  private void doRandomVehicleSwaps(int multiplier) {
    // randomly selects vehicles and trys to do an improving customer swap
    for (int k = 0; k < multiplier * V; k++) {
      int n = rng.nextInt(N); 
      int m = rng.nextInt(N);  
      if (n == 0 || m == 0 || n == m) continue;
      doRoutesSwap(n, m);
    }
  }
  
  private boolean trySwap(int n, int cn, int m, int cm) {
    // swap customer cn in vehicle n with customer cm in vehicle m if the total route gets shorter
    
    //System.out.println("try swap n, cn, m, cm = " + n + ", " + cn + ", " + m + ", " + cm);

    boolean swapped = false;
    Vehicle vn = vehicles[n];
    Vehicle vm = vehicles[m];
    
    // check capacities
    if (vn.spare() + customers[cn].demand() - customers[cm].demand() < 0) return false;
    if (vm.spare() + customers[cm].demand() - customers[cn].demand() < 0) return false;
    
    // set up vehicle n
    int Nvn = vn.numCustomers();
    double[] xn = new double[Nvn];
    double[] yn = new double[Nvn]; 
    int[] vnMap = new int[Nvn];
    for (int j = 0; j < Nvn; j++) {
      int k = (vn.route(j) == cn) ? cm : vn.route(j);
      xn[j] = customers[k].location().x();
      yn[j] = customers[k].location().y();
      vnMap[j] = k;
    }
    TspSolver tspn = new TspSolver(xn, yn, false, true, true, tspAnneal);
    
    // set up vehicle m    
    int Nvm = vm.numCustomers();
    double[] xm = new double[Nvm];
    double[] ym = new double[Nvm]; 
    int[] vmMap = new int[Nvm];
    for (int j = 0; j < Nvm; j++) {
      int k = (vm.route(j) == cm) ? cn : vm.route(j);
      xm[j] = customers[k].location().x();
      ym[j] = customers[k].location().y();
      vmMap[j] = k;
    }
    TspSolver tspm = new TspSolver(xm, ym, false, true, true, tspAnneal);
    
    double delta = (tspn.length() + tspm.length()) - (vn.length() + vm.length());
    
    if (delta < 0 || anneal(delta)) {
      // make the swap
      //System.out.println("swapping");
      
      swapped = true;
      
      Customer c = customers[cm];
      vm.decrement(c.demand());
      vn.increment(c.demand());
      c.clearVehicle();
      c.setVindex(n);
      vn.setRoute(tspn.mappedTour(vnMap));
      
      c = customers[cn];
      vn.decrement(c.demand());
      vm.increment(c.demand());
      c.clearVehicle();
      c.setVindex(m);
      vm.setRoute(tspm.mappedTour(vmMap));
      
      //printSolution(true);
    }
    
    return swapped;
  }
  
  private void doRandomChainReplacements(int multiplier) {
    // perform random 3-step chain replacements
    // a replacement gerneralizes a swap
    // instead of c1 <-> c2, c1 -> c2 -> c3 -> c1
    
    for (int k = 0; k < multiplier * N; k++) {
      
      int c1 = rng.nextInt(N);    // integer in the range [1,N)
      int c2 = rng.nextInt(N); 
      int c3 = rng.nextInt(N); 
      
      if (c1 == 0 || c2 == 0 || c3 == 0) continue;       // customers must not be 0 (warehose)
      if (c1 == c2 || c2 == c3 || c3 == c1) continue;   // customers must be distinct
      
      tryChainReplacement(c1, c2, c3);
      tryChainReplacement(c2, c1, c3);
      tryChainReplacement(c3, c2, c1);

    }
  }
  
  private void tryChainReplacement(int c1, int c2, int c3) {
    // c1, c2, c3 are distinct, non-zero cindexes
    
    int v1 = customers[c1].vindex();  // customer's current assigned vehicle
    int v2 = customers[c2].vindex();
    int v3 = customers[c3].vindex();
    
    if (v1 == v2 || v2 == v3 || v3 == v1) return;   // customers must be distinct
    
    // will the replacements fit?
    if (vehicles[v1].spare() + customers[c1].demand() - customers[c3].demand() < 0 ) return;
    if (vehicles[v2].spare() + customers[c2].demand() - customers[c1].demand() < 0 ) return;
    if (vehicles[v3].spare() + customers[c3].demand() - customers[c2].demand() < 0 ) return;
    
    // have 3 distinct customers currently assigned to distinct vehicles
    // set up the swap tsps
    
    
    // set up vehicle n
    TspSolver tsp1 = replaceTsp(v1, c1, c3);      
    TspSolver tsp2 = replaceTsp(v2, c2, c1);      
    TspSolver tsp3 = replaceTsp(v3, c3, c2); 
    
    double delta = vehicles[v1].length() + vehicles[v2].length() + vehicles[v3].length();
    double deltaTsp = tsp1.length() + tsp2.length() + tsp3.length();
    if (deltaTsp < delta || anneal(deltaTsp - delta)) {
      // make the replacements
      
      vehicles[v1].decrement(customers[c1].demand());
      vehicles[v2].decrement(customers[c2].demand());
      vehicles[v3].decrement(customers[c3].demand());
      vehicles[v1].increment(customers[c3].demand());
      vehicles[v2].increment(customers[c1].demand());
      vehicles[v3].increment(customers[c2].demand());  
      customers[c1].clearVehicle();
      customers[c2].clearVehicle();      
      customers[c3].clearVehicle();      
      customers[c1].setVindex(v2);
      customers[c2].setVindex(v3);      
      customers[c3].setVindex(v1);
      vehicles[v1].setRoute(tsp1.mappedTour());     
      vehicles[v2].setRoute(tsp2.mappedTour());     
      vehicles[v3].setRoute(tsp3.mappedTour()); 
            
    }    
  }
  
  private TspSolver replaceTsp(int vindex, int c1, int c2) {
    // returns a TspSolver object associated with replacing customer c1 in vehicle v with customer c2
    
    Vehicle v = vehicles[vindex];
    int Nv = v.numCustomers();
    double[] x = new double[Nv];
    double[] y = new double[Nv]; 
    int[] vMap = new int[Nv];
    for (int j = 0; j < Nv; j++) {
      int k = (v.route(j) == c1) ? c2 : v.route(j);
      x[j] = customers[k].location().x();
      y[j] = customers[k].location().y();
      vMap[j] = k;
    }
    TspSolver tsp = new TspSolver(x, y, false, true, true, tspAnneal);  
    tsp.setMap(vMap);
    return tsp;  
  }
    
  private void nnGreedyRoute(int vindex, int[] route, Bag<Integer> customerBag) {
    // greedy solution for vehicle vindex: visit next closest unvisited customer (nearest neighbor)
    // until full or no more customers
    
    Vehicle v = vehicles[vindex];
    int numCustomers = 1;
    int first = 0;
    route[0] = first;
    int prev = first;
    for (int i = 1; i <= route.length; i++) {
      int next = -1;
      double mind = Double.MAX_VALUE;
      for (int j : customerBag) {
        if (visited[j]) continue;                   // skip already visited customers
        if (customers[j].demand() <= v.spare()) {    // only look for customers that fit
          Point p = customers[j].location();
          Point q = customers[prev].location();
          
          // is point j outside the outer enclosing box for prev, ie not closer than min point?
          if (p.x() > q.x() + mind || p.x() < q.x() - mind) continue;
          if (p.y() > q.y() + mind || p.y() < q.y() - mind) continue;
                    
          double d = p.dist(q);
          if (d < mind) {
            mind = d;
            next = j;
          }
        }
      }

      if (next == -1) break;
      
      route[i] = next;
      v.increment(customers[next].demand());
      customers[next].setVindex(vindex);
      visited[next] = true;
      numCustomers++; 
      numVisited++;
      prev = next;
    }

    if (numCustomers == 0) {
      v.setInactive();
    } else {
 
      v.setNumCustomers(numCustomers);
      v.setRoute(route);
      
    }
    
  }
    
  private boolean addUnassignedToRoute() {
    int cindex = -1;
    for (int j = 1; j < N; j++) {
      if (!customers[j].assigned()) {
        cindex = j;
        break;
      }
    }
    if (cindex == -1) return false;
    
    Customer c = customers[cindex];
    for (int i = 0; i < V; i++) {
      // check the vehicles
      Vehicle v = vehicles[i];
      //System.out.println("i, demand, spare = " + i + " " + c.demand() + " " + v.spare());
      if (v.active()) {
        if (c.demand() <= v.spare()) {
          // v can handle c, add it
          v.add(c.cindex());
          numVisited++;
          return true;
        } else {
          // does v have an assigned customer with *less* demand than c?
          int[] route = v.route();
          for (int j = 1; j < v.numCustomers(); j++) {
            Customer crj = customers[route[j]];
            if (crj.demand() < c.demand()) {
              if (v.spare() + crj.demand() >= c.demand()) {
                //System.out.println("swap add = " + v.vindex() + ": " + c.cindex() + " for " + crj.cindex());
                v.remove(crj.cindex());
                //crj.clearVehicle();
                v.add(c.cindex());
                return true;
              }
              
            }
          }
          
        }
      }
    }
    return false;
  }
  
  private void addUnassignedToRoute(Bag<Integer> bag) {
    // a recursive version that "pushes" customers into adjacent vehicles
    for (int j : bag) {
      if (j == 0) continue;
      assignToNextVehicle(j, 0);
    }
  }
  
  private boolean assignToNextVehicle(int j, int i) {
    // assign customer j to vehicle i 
    // if necessary push a customer currently assigned to i to the next vehicle
    
    if (i == V) i = 0;
    
    Customer cj = customers[j];
    Vehicle vi = vehicles[i];
    
    if (cj.demand() <= vi.spare()) {
      // j fits, add it to i
      vi.add(j);
      numVisited++;
      return true;
      
    } else {
      // j does not fit, need to push some current customer in i to i+1
      int[] route = vi.route();
      for (int ji : route) {
        if (ji == 0) continue;    // skip the "warehouse" customer
        if (vi.spare() + customers[ji].demand() - cj.demand() >= 0) {
          vi.remove(ji);
          vi.add(j);
          return (assignToNextVehicle(ji, i+1));
        }
      }
    }
    
    return false;
    
  }
  
  private void addToNearestRadialRoute(int j) {
    // add customer j to nearest radial route regardless of vehicle capacity
    // nearest radial route is determined by the route of j's closest radial neighbor
    
    double thetaj = (new Theta(j)).theta();
    
    double min = Double.MAX_VALUE;
    int minj1 = -1;
    for (int j1 = 1; j1 < N; j1++) {
      if (!customers[j1].assigned()) continue;
      double thetaj1 = (new Theta(j1)).theta();
      double delta = Math.abs(thetaj - thetaj1);
      if (delta < min) {
        min = delta;
        minj1 = customers[j1].vindex();
      }
    }
    vehicles[minj1].add(j);
    numVisited++;
  }
  
  
  
  /*
  *   clustering
  */
  
  private int calcK() {
    int k = 0;                              // number of clusters - todo: get clever about K
    double totalDemand = 0;
    for (int j = 0; j < N; j++) totalDemand += customers[j].demand();
    k = (int) (totalDemand / vehicles[0].capacity()) + 1;
    if (k > V) k = V; 
    return k;   
  }
  
  private int[] kMeansClusters(int K) {
    // returns an array of bag objects, each bag representing the customers in one of K clusters
    // uses k-means clustering
    
    // compute the initial centroids
    // 1. x-y ranges
    double xbar = 0;
    double ybar = 0;
    double maxx = Double.MIN_VALUE;
    double minx = Double.MAX_VALUE;
    double maxy = Double.MIN_VALUE;
    double miny = Double.MAX_VALUE; 
    for (int j = 0; j < N; j++) {
      Point p = customers[j].location();
      xbar += p.x();
      ybar += p.y();
      if (p.x() < minx) minx = p.x();
      if (p.x() > maxx) maxx = p.x();
      if (p.y() < miny) miny = p.y();
      if (p.y() > maxy) maxy = p.y();
    }  
  
    // 2. random dispersal of initial centroids in the x-y ranges
    Point[] centroids = new Point[K];

    for (int k = 0; k < K; k++) {
      double x = minx + rng.nextDouble() * (maxx - minx);
      double y = miny + rng.nextDouble() * (maxy - miny);
      centroids[k] = new Point(x, y);      
    }

    return kMeansClusters(centroids);
    
  }
  
  private int[] kMeansClusters(Point[] centroids) {
    // converge the cluster centroids
    K = centroids.length;
    int[] clusters = new int[N];    // cluster[j] is the cluster containing customer j
    
    boolean looping = true;
    int maxLoops = 100;
    int numLoops = 0;
    while (looping && numLoops < maxLoops) {
      // assign each customer to its nearest centroid/cluster
      for (int j = 1; j < N; j++) {     // customer 0 is the warehouse, don't cluster it
        Customer c = customers[j];
        double min = Double.MAX_VALUE;
        for (int k = 0; k < K; k++) {
          double d = centroids[k].dist(c.location());
          if (d < min) {
            min = d;
            clusters[j] = k;
          }
        }
      }
      
      // recompute each cluster centroid
      double[] sumx = new double[K];
      double[] sumy = new double[K];
      int[] count = new int[K];
      double delta = 0;

      for (int j = 0; j < N; j++) {
        sumx[clusters[j]] += customers[j].location().x();
        sumy[clusters[j]] += customers[j].location().y();
        count[clusters[j]]++;
      }
      for (int k = 0; k < K; k++) {
        sumx[k] /= count[k];
        sumy[k] /= count[k];
        Point p = new Point(sumx[k], sumy[k]);
        delta += centroids[k].dist(p);
        centroids[k] = p; 
               
      }
      
      //converged?
      double eps = .001;
      if (delta < eps) looping = false;
      numLoops++;
    }
    
    return clusters;    
  }
  

  private int[] angleSweepClusters() {    
    // build a minPQ of thetas (see Theta object below)
    // note: customer 0 is the warehouse and not in thetasPQ
    // "sweep" through the customers in min theta order, adding to a vehicle as long as capacity
    
    int[] clusters = new int[N];
 
    MinPQ<Theta> thetasPQ = new MinPQ<Theta>();
    for (int j = 1; j < N; j++) thetasPQ.insert(new Theta(j));
    int vindex = 0;
    double spare = capacity;
    while (!thetasPQ.isEmpty()) {
      int cindex = thetasPQ.delMin().cindex();
      if (customers[cindex].demand() > spare) {
        ++vindex;
        spare = capacity;
      }
      spare -= customers[cindex].demand();
      clusters[cindex] = vindex;
    }
    K = (vindex + 1 > V) ? V : vindex + 1;

    return clusters;
    
  }
  
  private int[] angleSweepClusters(double thetaOffset) {
    // build a minPQ of thetas (see Theta object below)
    // note: customer 0 is the warehouse and not in thetasPQ
    // "sweep" through the customers in min theta order, adding to a vehicle 
    // number of vehicles is given as k

    //System.out.println(thetaOffset);
    
    int[] clusters = new int[N]; 
 
    double twoPi = 2 * Math.PI;
    MinPQ<Theta> thetasPQ = new MinPQ<Theta>();
    for (int j = 1; j < N; j++) {
      Theta t = new Theta(j);
      if (t.theta() > twoPi - thetaOffset) t.setTheta(-(twoPi-t.theta()));
      //System.out.println("customer, point, theta = " + customers[j].cindex() + ", " + customers[j].location() + ", " + t.theta());
      thetasPQ.insert(t);
    }
    int vindex = 0;
    double spare = capacity;
    while (!thetasPQ.isEmpty()) {
      Theta t = thetasPQ.delMin();
      int cindex = t.cindex();
      //System.out.println("cindex, theta = " + cindex + ", " + t.theta());
      if (customers[cindex].demand() > spare) {
        ++vindex;
        spare = capacity;
      }
      //System.out.println("customer, location, theta, cluster = " + cindex + ", " + customers[cindex].location() + ", " + t.theta() + ", " + vindex);
      spare -= customers[cindex].demand();
      clusters[cindex] = vindex;
    }
    K = (vindex + 1 > V) ? V : vindex + 1;

    return clusters;
  }
  
  /*
  * some misc and helper methods
  */
  
  private boolean feasibleRoute(boolean verbose) {
    // is the total routing feasible?
    boolean[] served = new boolean[N];
    boolean feasible = true;
    for (int i = 0; i < V; i++) {
      Vehicle v = vehicles[i];
      if (v.spare() < 0) {
        if (verbose) System.out.println("vehicle over capacity (spare): " + v.vindex() + " (" + v.spare() + ")");
        feasible = false;
      }
      for (int j = 0; j < v.numCustomers(); j++) served[customers[v.route(j)].cindex()] = true;
    }
    for (int j = 0; j < N; j++) {
      if (!served[j]) {
        if (verbose) System.out.println("customer not served = " + j);
        feasible = false;
      }
    }
    return feasible;
  }
  
  
  private double printSolution(boolean verbose) {
    // print the solution in the specified output format        
    double length = totalRoutesLength();
    System.out.println(length + " 0");
    if (verbose) {
      for (int i = 0; i < V; i++) System.out.println(vehicles[i]);
    }

    return length;
  }
  
  private void printForPlot() {
    feasibleRoute(false);
    System.out.println("\n" + totalRoutesLength());
    String s = "";
    for (int i = 0; i < V; i++) s += vehicles[i] + " ";
    System.out.println(s + "\n");
  }
  
  private double totalRoutesLength() {
    double length = 0;
    for (int i = 0; i < V; i++) length += vehicles[i].length();
    return length;
  }
  
  private double dist(int i, int j) {
    // helper method: distance from customers[i] to customers[j]
    return customers[i].distTo(customers[j]);
  }
  
  /*
  *  the vehicle model
  */

  private class Vehicle {
    private int vindex;
    private double capacity, used;
    private boolean active;
    private int[] route;
    private int numCustomers;
    private Vehicle (int i, double c) {
      vindex = i;
      capacity = c;
      active = false;
      numCustomers = 0;
      route = new int[N];
    }
    private int vindex() { return vindex; }
    private double capacity() { return capacity; }
    private boolean active() { return active; }
    private void setActive() { 
      active = true; 
      numCustomers = 0;
    }
    private void setInactive() { 
      active = false;
      used = 0;
      for (int j = 1; j < numCustomers; j++) customers[route[j]].clearVehicle();
      numCustomers = 0; 
      route = new int[N];
    }
    private double spare() { return capacity - used; }
    private void increment(double d) { used += d; }
    private void decrement(double d) { used -= d; }
    private int[] route() { return route; }
    private int route(int j) { return route[j]; }
    private void setRoute(int[] r) { route = r; }
    private int numCustomers() { return numCustomers; }
    private void setNumCustomers(int s) { numCustomers = s; }
    private double length() {
      if (!active) return 0;
      double length = 0;
      for (int j = 0; j < numCustomers-1; j++) {
        length += customers[route[j]].location().dist(customers[route[j+1]].location());
      }
      length += customers[route[numCustomers-1]].location().dist(customers[route[0]].location());
      return length;
    }
    private boolean remove(int cindex) {
      // removes customer from vehicles route
      // todo: optimize the route after removal
      boolean cindexFound = false;
      for (int j = 1; j < numCustomers - 1; j++) {
        if (cindexFound || route[j] == cindex) {
          route[j] = route[j+1];
          cindexFound = true;
        }
      }
      if (cindexFound) {
        customers[cindex].clearVehicle();
        decrement(customers[cindex].demand());
        numCustomers--;
        //System.out.println("removed customer from vehicle = " + cindex + ", " + vindex);
        
      }
      return cindexFound;
    }
    private boolean add(int cindex) {
      // adds customer to vehicle's route
      // todo: optimize the route after addition
      for (int j = 1; j < numCustomers; j++) {
        if (route[j] == cindex) return false;      // already in route, can't add again
      }
      if (++numCustomers > route.length) {
        //need to expand the route array size
        int[] temp = new int[numCustomers];
        for (int j = 0; j < numCustomers - 1; j++) temp[j] = route[j];
        route = temp;
      }
      route[numCustomers-1] = cindex;
      increment(customers[cindex].demand());
      customers[cindex].setVindex(this.vindex);
      //System.out.println("added customer to vehicle = " + cindex + ", " + vindex);
      return true;
    }
    private boolean insert(int cindex) {
      // inserts a customer into an existing vehicle route so that total route length is minimized
      
      // degenerate case
      if (numCustomers == 1) {
        return add(cindex);
      } 
      
      // is customer already on the route?
      for (int j = 1; j < numCustomers; j++) {
        if (route[j] == cindex) return false;      // already in route, can't add again
      }  
      
      // where should customer go?
      double l = length();
      double min = Double.MAX_VALUE;
      int minj = -1;
      for (int j = 0; j < numCustomers; j++) {
        int jp1 = ((j + 1) < numCustomers) ? j+1 : 0;
        double test = l - dist(route[j], route[jp1]) + dist(route[j], cindex) + dist(cindex, route[jp1]);
        if (test < min) {
          min = test;
          minj = j;
        }
      } 
      
      if (minj == -1) return false;   // unexpected 
      
      // insert the customer
      for (int j = minj+1; j < numCustomers; j++) {
        route[j+1] = route[j];
      } 
      route[minj + 1] = cindex;
      numCustomers++;
      
      return true;
      
    }
    public String toString() {
      // returns string represenaation of the route formatted for output
      if (!active) return "0 0";
      String s = ""; 
      for (int i = 0; i < numCustomers; i++) { s += route[i] + " "; }
      return s + "0";
    }
    
  }
  
  
  /*
  *  the customer model
  */
  
  private class Customer implements Comparable<Customer>{
    private int cindex, vindex;
    private double demand;
    private Point location;
    private Customer (int j, double d, Point p) {
      cindex = j;
      vindex = -1;
      demand = d;
      location = p;
    }
    private int cindex() { return cindex; }
    private int vindex() { return vindex; }
    private void setVindex(int i) { vindex = i; }
    private double demand() { return demand; }
    private Point location() { return location; }
    private boolean assigned() { return vindex != -1; }
    private void clearVehicle() { vindex = -1; }
    public int compareTo(Customer that) {
      if (demand > that.demand()) return 1;
      if (demand < that.demand()) return -1;
      return 0;
    }
    private double distTo(Customer that) {
      // distance to the given customer
      return this.location().dist(that.location());
    }
    public String toString() {
      return "customer, demand : " + cindex + ", " + demand;
    }
  }

  private class Point {
    private double x, y;
    private Point (double a, double b) {
      x = a;
      y = b;
    }
    private double x() { return x; }
    private double y() { return y; }
    private double dist(Point that) {
      // distance from this point to that point
      double tx = x - that.x();
      double ty = y - that.y();
      return Math.sqrt(tx * tx + ty * ty);
    }
    public String toString() {
      return "(" + x + ", " + y + ")";
    }
  }
  
  private class Theta implements Comparable<Theta> {
    private double theta;
    private int cindex;
    private Theta (int c) {
      cindex = c;
      Point w = customers[0].location();
      Point p = customers[cindex].location();
      theta = Math.atan2(p.y() - w.y(), p.x() - w.x());
      if (theta < 0) theta = 2 * Math.PI + theta;
      //System.out.println("cindex, theta = " + cindex + ", " + theta);
    }
    private double theta() { return theta; }
    private void setTheta(double t) { theta = t; }
    private int cindex() { return cindex; }
    public int compareTo (Theta that) {
      if (theta > that.theta()) return 1;
      if (theta < that.theta()) return -1;
      return 0;
    }
  }

    
}