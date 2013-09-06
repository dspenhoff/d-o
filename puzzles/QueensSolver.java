import java.util.Random;

class QueensSolver {

  private int N;
  private int[] q;          // q[j] is the row for the queen in queen j
  private int[] conflicts;  // conflicts[j] is the number of conflicts for queen j
  
  private QueensSolver(int n) {
    

    N = n;
    q = new int[N];
    if (!solverB()) {
      System.out.println("no solution found.");
      printSolution();
    } else {
      printSolution();
    }
  }
  
  /*
  * solverA is exhaustive search using recursion
  * uses the setQueen and blocked methods
  */
  
  private boolean solverA() {
    return setQueen(0);
  }
  
  private boolean setQueen(int col) {
    // place queen in column given placements in columns j < col
    if (col >= N) return true;
    for (int row = 0; row < N; row++) {
      if (!blocked(row, col)) {
        q[col] = row;
        if (setQueen(col + 1)) return true;
      }
    }
    return false;
  }
  
  private boolean blocked(int row, int col) {
    // is square (row, col) blocked given q[j] for j < col
    for (int j = 0; j < col; j++) {
      if (q[j] == row || q[j] == row - (col - j) || q[j] == row + (col - j)) return true;
    }
    return false;
  }
  
  /*
  * solverB is a local search to eliminate conflicts in an initial solution
  */
  
  private boolean solverB() {
    // 1. randomly assign queens
    Random rng = new Random();
    for (int j = 0; j < N; j++) q[j] = rng.nextInt(N);
    
    // 2. count the conflicts
    conflicts = new int[N];
    for (int j = 0; j < N; j++) conflicts[j] = numConflicts(q[j], j);
    
    // 3. local search loop:
    //    pick the queen with the most conflicts
    //    move it to it's minimum conflict position
    //    update the conflicts
    boolean done = false;
    int doneCount = 0;
    while (!done && doneCount < 10000) {
      int max = Integer.MIN_VALUE;
      int maxq = 0;
      for (int j = 1; j < N; j++) {
        if (conflicts[j] > max) {
          max = conflicts[j];
          maxq = j;
        }
      }
      
      int min = Integer.MAX_VALUE;
      int newq = 0;
      for (int i = 0; i < N; i++) {
        if (i == q[maxq]) continue;
        int n = numConflicts(i, maxq);
        if (n <= min) {
          min = n;
          newq = i;
        }
      }
      q[maxq] = newq;
      int sum = 0;
      for (int j = 0; j < N; j++) {
        int n = countConflicts(j);
        sum += n;
        conflicts[j] = n;
      }
      
      printQueens();
      printConflicts();
      if (sum == 0) done = true;
      doneCount++;
    }
    

    
    return false;
  }
  
  private int countConflicts(int col) {
    // number of conflicts for queen in column
    return numConflicts(q[col], col); 
  }
  
  private int numConflicts(int row, int col) {
    int count = 0;
    for (int j = 0; j < N; j++) {
      if (j == col) continue;
      if (q[j] == row)                  count++; 
      else if (q[j] == row - (col - j)) count++; 
      else if (q[j] == row + (col - j)) count++; 
    } 
    return count;   
  }
  
  private void printQueens() {
    String s = "queens: ";
    for (int j = 0; j < N; j++) s += q[j] + " ";
    System.out.println(s);
  }
  
  private void printConflicts() {
    String s = "conflicts: ";
    for (int j = 0; j < N; j++) s += conflicts[j] + " ";
    System.out.println(s);
  }
  
  private void printSolution() {
    System.out.println(N);
    String s = "";
    for (int j = 0; j < N; j++) s += q[j] + " ";
    System.out.println(s);
  }
  
  public static void main(String[] args) {
    //System.out.println("QueensSolver " + args[0]);
    new QueensSolver(Integer.parseInt(args[0]));
  }
}