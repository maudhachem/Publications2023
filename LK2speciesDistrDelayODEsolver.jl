# Filename: LK2speciesDistrDelayODEsover.jl
# Author: Maud El-Hachem and Nick Beeton, CSIRO Hobart, March 2023
# Description: This file contains the code to solve equations (22)-(27)
# and displays the solutions N1(t) and N2(t) shown in Figure 2(a)-(c)
# from the following reference
# M El-Hachem and NJ Beeton (2023) 
# Coexistence in two-species competition with delayed maturation

using DifferentialEquations
using Plots
using LaTeXStrings
using CSV
using DataFrames
using Random, Distributions
# This function solves equations (22)-(27)
# Input: 	du is the vector at the LHS of equations (22)-(27)
#			u is the vector containing N1, N2, f1, f2 x1_1 to x1_p1, x2_2 to x2_p2,
#			y1_1 to y1_p1 and y2_2 to y2_p2 (22)-(27)
#			param are the parameters used are given in the following order
#			m1, m2, r1, r2, mu1, mu2, a1, a2, b1, b2, beta1, beta2, p1, p2, Delta1, Delta2,
#			alpha1, alpha2
#			t is the vector containing the beginning and the ending time
function solve_odesys(du,u,param,t)
	#parameters
	m1, m2, r1, r2, mu1, mu2, a1, a2, b1, b2, beta1, beta2, p1, p2, Delta1, Delta2, alpha1, alpha2 = param
	
	#N1(t) = u[1]
	#N2(t) = u[2]
	#f1(t) = u[3]
	#f2(t) = u[4]
	du[1] = r1*(beta1/(mu1+beta1))^p1*u[p1-1+5]*exp(-Delta1*(a1*u[3] + b2 * u[4]))-u[1]*m1
	du[2] = r2*(beta2/(mu2+beta2))^p2*u[2*(p1-1)+7+p2-1]*exp(-Delta2*(a2 * u[4] + b1 * u[3]))-u[2]*m2
	du[3] = (u[1] - u[2*(p1-1)+6])/Delta1
	du[4] = (u[2] - u[2*(p1-1)+8+2*(p2-1)])/Delta2
	#x1_1 = u[5]
	du[5] = (beta1+mu1)/alpha1*(u[1]-u[5])
	#x1_2 to x1_p1
	for i in 2:p1 
		j = i-1
		du[j+5] = (beta1+mu1)/alpha1*(u[j+5-1]-u[j+5]) 
	end
	#y1_1
	du[p1-1+6] = beta1/alpha1*(u[1]-u[p1-1+6])
	#y1_2 to y1_p1
	for i in 2:p1
		j = i - 1
		du[p1-1+6+j] = (beta1)/alpha1*(u[p1-1+6+j-1]-u[p1-1+6+j])
	end
	#x2_1
	du[2*(p1-1)+7] = (beta2+mu2)/alpha2*(u[2]-u[2*(p1-1)+7])
	#x2_2 to x2_p1
	for i in 2:p2
		j = i - 1
		du[2*(p1-1)+7+j] = (beta2+mu2)/alpha2*(u[2*(p1-1)+7+j-1]-u[2*(p1-1)+7+j])
	end
	#y2_1
	du[2*(p1-1)+8+p2-1] =(beta2)/alpha2*(u[2]-u[2*(p1-1)+8+p2-1])
	#y2_2 to y2_p2
	for i in 2:p2
		j = i - 1
		du[2*(p1-1)+8+p2-1+j] = (beta2)/alpha2*(u[2*(p1-1)+8+p2-1+j-1]-u[2*(p1-1)+8+p2-1+j])
	end
end

# This function computes and displays the numerical solutions N1(t) and N2(t) shown in Figure 2(a)-(c)
# Input: 	parameters used equations (22)-(27) given in the following order
#			m1, m2, r1, r2, mu1, mu2, a1, a2, b1, b2, ptab, betatab
#			ptab is a 2x2 table: 1st line contains values of p1 to be tested,
#			2nd lines contains values of p2
#			betatab is a 2x2 table of values of beta1 or beta2 corresponding 
#			to each value of p1 or p2 from table ptab
# 			qt: the quantile of the original distribution indicating where the gamma distribution is truncated
#			strOutputFile: name of the svg file where to save the plot
function solve_display(m1, m2, r1, r2, mu1, mu2, a1, a2, b1, b2, ptab, betatab, qt, strOutputFile)
	
	tspan = (0.0,60.0)
	
	# size of the table of the shapes of the gamma distribution 
	tabsize = length(ptab)
	for index = 1:tabsize
		# shape of the gamma distribution for species 1 and 2	
		p1 = ptab[1][index]
		p2 = ptab[2][index]
		# rate of the gamma distribution for species 1 and 2
		beta1 = betatab[1][index]
		beta2 = betatab[2][index]

		# location of the truncation for species 1 and 2
		Delta1 = quantile(Gamma(p1,1/beta1), qt)
		Delta2 = quantile(Gamma(p2,1/beta2), qt)
		println(Delta1)
		println(Delta2)
		
		# total number of equations
		sizeu = 2*p1+2*p2+4
		
		#calculating the parameter alpha1 and alpha2
		if (p1==1)
			alpha1 = 1-exp(-beta1*Delta1)
		elseif (p1==2)
			alpha1 = 1-exp(-beta1*Delta1)*(1+(beta1*Delta1))
		elseif (p1>2)
			alpha1 = 1-exp(-beta1*Delta1)*sum(((beta1*Delta1)^n)/factorial(n) for n in 0:p1-1)
		end 
		if (p2==1)
			alpha2 = 1-exp(-beta2*Delta2)
		elseif (p2==2)
			alpha2 = 1-exp(-beta2*Delta2)*(1+(beta2*Delta2))
		elseif (p2>2)
			alpha2 = 1-exp(-beta2*Delta2)*sum(((beta2*Delta2)^n)/factorial(n) for n in 0:p2-1)
		end
#		println(alpha1)
#		println(alpha2)

		# initialisation of functions N1, N2, f1, f2, xk, yk for k=1..i and species 1 and 2
		u0 = zeros(sizeu)
		# N1(0)
		u0[1] = 0.5
		# N2(0)
		u0[2] = 0.5
		# f1(0)
		u0[3] = 0
		# f2(0)
		u0[4] = 0
		# xk(0), yk(0) for k=1..i and species 1 and 2
		for i in 5:sizeu
			u0[i] = 0
		end
		# parameters in the equations	
		param = (m1, m2, r1, r2, mu1, mu2, a1, a2, b1, b2, beta1, beta2, p1, p2, Delta1, Delta2, alpha1, alpha2)

		# define the ODE problem
		prob = ODEProblem(solve_odesys,u0,tspan,param)
		# call the solver BS3 to get the ODE solutions
		sol1 = solve(prob, BS3(), abstol = 1e-6, reltol = 1e-6)
		
		# ticks labels to display on the axis
		ticksy = ([0, 0.2, 0.4, 0.6, 0.8, 1.0],[L"0", L"0.2", L"0.4", L"0.6", L"0.8", L"1.0"])
		ticksx = ([0,20,40,60],[L"0",L"20",L"40",L"60"])

		# display the solutions
		if (index == 1)
		# display the solutions for p1=p2=1
			plot(sol1,idxs=[(0,1),(0,2)],ylim=(0, 1),xlim=(0,60),size = (350, 300),
			linewidth=2,linecolor = [:dodgerblue2 :red2 ])
		else
		# display the solutions for p1=p2=20 on the same plot as for p1=p2=1
			plot!(sol1,idxs=[(0,1),(0,2)],ylim=(0, 1),xlim=(0,60),size = (350, 300),
			labelfontsize=14,tickfontsize = 14,
			linewidth=2,linecolor = [:limegreen :goldenrod2],
			xticks = ticksx, 
			yticks=ticksy, 
			legend=false, 
			grid=false, 
			framestyle = :box,
			xlabel = L"t",
			ylabel = L"N_1(t) \quad N_2(t)")
		end
		# save the figure 
		strsvg = strOutputFile * ".svg"
		savefig(strsvg)
		
# uncomment the next three lines if you want to output the numerical solutions in a cvs file
#		df = DataFrame(sol1)
#		strcsv = strOutputFile * string(index) * ".csv"
#		CSV.write(strcsv,string.(df))
	end
end

# find and display the solutions of Figure 2(a)
# m1=m2=mu1=mu2=0.3, r1=r2=1, a1=a2=1, b1=b2=0.1, p1=1 and 20, p2=1 and 20, beta1 = 1 and 20, beta2 = 1/8 and 20/8, quantile = 0.99
solve_display(0.3, 0.3, 1, 1, 0.3, 0.3, 1, 1, 0.1, 0.1, [[1,20],[1,20]], [[1,20],[1/8,20/8]], 0.99, "Figure3a")
# find and display the solutions of Figure 2(b)
# m1=m2=mu1=mu2=0.3, r1=r2=1, a1=a2=1, b1=b2=0.1, p1=1 and 20, p2=1 and 20, beta1 = 1 and 20, beta2 = 1/2 and 10, quantile = 0.99
solve_display(0.3, 0.3, 1, 1, 0.3, 0.3, 1, 1, 0.1, 0.1, [[1,20],[1,20]], [[1,20,],[1/2,10]],  0.99, "Figure3b")
# find and display the solutions of Figure 2(c)
# m1=m2=mu1=mu2=0.3, r1=r2=1, a1=a2=1, b1=b2=3, p1=1 and 20, p2=1 and 20, beta1 = 1 and 20, beta2 = 1/2 and 10, quantile = 0.99
solve_display(0.3, 0.3, 1, 1, 0.3, 0.3, 1, 1, 3,   3,   [[1,20],[1,20]], [[1,20,],[1/2,10]],  0.99, "Figure3c")
