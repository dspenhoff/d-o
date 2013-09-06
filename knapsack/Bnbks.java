// Copyright (c) by David M. Spenhoff. All rights reserved

class Bnbks


	// implements a branch and bound algorithm to solve the 0-1 knapsack problem
	// Horowitz-Sahni algorithm taken from Martello and Toth, pp 30-31
	
 	public Bnbks(n, p, w, c)
 	/*
 		# following the terminology in M&T, p is the 'profit' array, w is the 'weight' array
 		# and c is the total capacity of the knapsack. p and w are asumed to be ordered such
 		# that p[j]/w[j] decreases as j increases. note that indexes are offset -1 (e.g., 0..(n-1) 
 		# instead of 1..n)
 		
 		# test data from M&T ... solution is x=1001000, z=107
 		# p = [70, 20, 39, 37, 7, 5, 10]
 		# w = [31, 10, 20, 19, 4, 3, 6]
 		# c = 50
 	*/
 		
 		double start_time = Time.now.to_f	# start timer
		
 		// 1. initialize
 		double z = z_hat = 0
 		c_hat = c
 		p[n] = 0
 		w[n] = 9999999
 		int x[] = new Array();
 		x_hat = Array.new
 		n.times do |k| x[k] = x_hat[k] = 0 end
 		j = 0

 		while (true) {	
 			// 2. compute upper bound u
 			sum_wk = sum_pk = 0
 			r = j
 			for (int k = j; k < n; k++) {
 				sum_wk += w[k]
 				sum_pk += p[k]
				if (sum_wk > c_hat) { 
					r = k
					break
				}		
 			}
 			u = (sum_pk - p[r]) + ((c_hat - (sum_wk - w[r])) * p[r].to_f / w[r]).floor

			if (z >= z_hat + u) {
				// no better solution is on this branch
				backtrack=true
			} else {
 				while (true)
	 				# 3. perform forward step
	 				while w[j] <= c_hat
						c_hat -= w[j]
						z_hat += p[j]
						x_hat[j] = 1
						j += 1
	 				end

	 				if j <= (n - 1)
	 					x_hat[j] = 0
	 					j += 1
	 				end
	 				
	 				# branch to next step
	 				if j < (n - 1)
	 					backtrack = false
	 					break		# go to 2
	 				elsif j == (n - 1)
	 					backtrack = false
	 					next		# go to 3
	 				else
	 					# 4. update the best solution so far
						if z_hat > z
							z = z_hat
							0.upto(n - 1) do |k|
								x[k] = x_hat[k]
							end
						end
						j = n - 1
						if x_hat[j] == 1
							c_hat += w[j]
							z_hat -= p[j]
							x_hat[j] = 0
						end
						backtrack = true
						break
					end
				end
			end
			
			if backtrack 
				# 5. backtrack
				i = j - 1
				while i >= 0
					if x_hat[i] == 1 then break end
					i -= 1
				end
				if i < 0 then break end		# exit - optimal solution found: x
				c_hat += w[i]
				z_hat -= p[i]
				x_hat[i] = 0
				j = i + 1
			end
		end

		# output test data to console		
		# puts "x = " + x.to_s
		# puts "z = " + z.to_s
		
		@elapsed_time += Time.now.to_f - start_time
		
		x				# return the optimal vector of selected items
	end
	
	# data accessor methods
	def elapsed_time
  	(@elapsed_time * 10000).to_i / 10.0
	end
	
end
