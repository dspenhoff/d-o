import java.util.Random;

class AllIntervalSeriesSolver {

  private int N;
  private int[] series;          
  private int[] deltas; 
  int[] deltaCounts; 
  
  private AllIntervalSeriesSolver(int n) {
    

    N = n;
    //series = new int[N];
    //deltas = new int[N];
    
    if (!solverB()) {
      System.out.println("no solution found.");
      printSolution();
    } else {
      printSolution();
    }
  }
  
  /*
  * solverA is exhaustive search using recursion
  */
  
  private boolean solverA() {
    return setSeries(0);
  }
  
  private boolean setSeries(int n) {
    // place queen in column given placements in columns j < col
    if (n >= N) return true;
    for (int k = 0; k < N; k++) {
      //System.out.println("n, k = " + n + ", " + k);
      if (n == 0) {
        series[n] = k;
        if (setSeries(n + 1)) return true;
      } else {
        if (!blocked(k, n)) {
          series[n] = k;
          deltas[n-1] = Math.abs(series[n-1] - series[n]);     
          if (setSeries(n + 1)) return true;
        }
      }
    }
    return false;
  }
  
  private boolean blocked(int k, int n) {
    // given a valid solution[0.. n-1], is solution[n] = k valid?
    int delta = Math.abs(k - series[n - 1]);
    for (int j = 0; j < n; j++) {
      if (series[j] == k || deltas[j] == delta) return true;
    }
    return false;
  }
  
  /*
  * solverB starts with a randomized series of size N and uses local search to move towards valid solution
  */
  
  private boolean solverB() {
    
    initializeSeries();
    
    Random rng = new Random();
    int numConflicts = conflicts();
    int maxSwaps = *;
    int numSwaps = 0;
    while (numSwaps < maxSwaps) {
      for (int i = 0; i < N; i++) {
        if (!hasConflicts(i)) continue;
        for (int j = 0; j < N; j++) {
          if (j == i || !hasConflicts(j)) continue;
          swap(i, j);
          numSwaps++;
          int conflicts = conflicts();
          if (conflicts >= numConflicts) {
            swap(j, i);
            numSwaps--;
          } else {
            numConflicts = conflicts;
            break;
          }
          if (valid()) return true;
          break;
        }
      }

    }
    
    return false;
    
  }
  
  private void swap(int i, int j) {
    // swaps series elements i and j and updates deltas
    int temp = series[i];
    series[i] = series[j];
    series[j] = temp;
    updateDelta(i);
    updateDelta(j);
  }
  
  private boolean hasConflicts(int i) {
    if (deltaCounts[deltas[i]] > 1) return true;
    return false;
  }
  
  private void updateDelta(int i) {
    if (i > 0) deltas[i-1] = Math.abs(series[i-1] - series[i]);
    if (i < N-1) deltas[i] = Math.abs(series[i] - series[i+1]);
    deltaCounts = new int[N];
    for (int j = 0; j < N-1; j++) deltaCounts[deltas[i]]++;
  }
  
  private int conflicts() {
    // returns the number of conflicts in the current series/deltas
    int conflicts = 0;
    deltaCounts = new int[N];
    for (int i = 0; i < N-1; i++) deltaCounts[deltas[i]]++;
    for (int j = 0; j < N; j++) if (deltaCounts[j] > 1) conflicts++;
    return conflicts;
  }
  
  private boolean valid() {
    int conflicts = 0;
    int[] counts = new int[N];
    for (int i = 0; i < N-1; i++) counts[deltas[i]]++;
    for (int j = 0; j < N; j++) if (counts[j] > 1) return false;
    return true;    
  }
  
  private void initializeSeries() {
    series = new int[N];
    deltas = new int[N-1];
    
    // randomize initial series
    for (int i = 0; i < N; i++) series[i] = i;
    Random rnd = new Random();
    for (int i = N-1; i >= 0; i--) {
      int k = rnd.nextInt(i + 1);
      int temp = series[k];
      series[k] = series[i];
      series[i] = temp;
    }

    // initialize deltas
    deltas = new int[N];
    for (int i = 0; i < N-1; i++) deltas[i] = Math.abs(series[i] - series[i+1]);
  }
 
  /*
  * output methods
  */
  
  private void printSolution() {
    System.out.println(N);
    printSeries();
  }
  
  private void printSeries() {
    String s = "";
    for (int j = 0; j < N; j++) s += series[j] + " ";
    System.out.println(s);    
  }
  
  private void printDeltas() {
    String s = "";
    for (int j = 0; j < N-1; j++) s += deltas[j] + " ";
    System.out.println(s);    
  }
  
  public static void main(String[] args) {
    new AllIntervalSeriesSolver(Integer.parseInt(args[0]));
  }
}