#===================================================================================
#  Lots of R packages for graphs
#	these are described in a CRAN task view here:
#	http://cran.r-project.org/web/views/gR.html
#	-we willl use the 'network' package
#===================================================================================

#install.packages("network")
library(network)
#help("network-package")

#===================================================================================
#	-Simulated population parameters
#===================================================================================

L<-num_obs<-200  #number of households and individuals
S<-rep(1,L)
E<-rep(0,L)
I<-rep(0,L)
R<-rep(0,L)


#========================================================================
#	Model two neighborhoods
#	-LD: general background connectivity
#	-LD2 & LD# within neighborhood extra connectivity
#========================================================================

connectivity<-.01

LD<-matrix(rbinom(L^2,1,connectivity),L)
for( i in 1:L)
    {
    LD[i,]<-LD[,i]
    }

diag(LD)<-0

connectivity<-.09

LD2<-matrix(rbinom((L/2)^2,1,connectivity),L/2)
for( i in 1:L/2)
    {
    LD2[i,]<-LD2[,i]
    }

diag(LD2)<-0

LD[1:(L/2),1:(L/2)]<-LD[1:(L/2),1:(L/2)]+LD2

connectivity<-.09
LD3<-matrix(rbinom((L/2)^2,1,connectivity),(L/2))
for( i in 1:L/2)
    {
    LD3[i,]<-LD3[,i]
    }

diag(LD3)<-0

LD[101:200,101:200]<-LD[101:200,101:200]+LD3

NET<-network(LD, directed=FALSE)
plot(NET)


#===================================================================================
#	Plotting functions
#	-COORDS lets it draw network once and saves the coordinates of each node
#	-COLS a vector to define colors in the plot 
#	-visualize: plots network, coloring the I nodes red
#===================================================================================


COORDS<-plot(NET,vertex.col="black")
COLS<-rep(1,L)

visualize_net<-function()
	{
	COLS[I==1]<-"red"
	COLS[E==1]<-"orange"
	COLS[S==1]<-"black"
	COLS[R==1]<-"white"
	plot(NET,vertex.col=COLS, jitter=FALSE, coord=COORDS)
	}


visualize_nodes<-function()
	{
	COLS[I==1]<-"red"
	COLS[E==1]<-"orange"
	COLS[S==1]<-"black"
	COLS[R==1]<-"white"
    lines(COORDS,col=COLS,pch=19,type="p") # much faster than points()
	}

#===================================================================================
#	-Set up Simulation 
#	-S,E,I,R vectors
#========================================================================

setup<-function()
	{
	S<<-rep(1,L)
	E<<-rep(0,L)
	I<<-rep(0,L)
	R<<-rep(0,L)
	recoverday<<-rep(0,L)
	}

#========================================================================
#	Key functions of stochastic simulation
#	-infect changes I from 0 to 1, also changes E and S to 0
#========================================================================

infect<-function(node,day=i)
	{
	I[node]<<-1
	E[node]<<-0
	S[node]<<-0
	cases<<-which(I==1)
	recoverday[node]<<-day+duration[node]
	}
		
expose<-function()
	{
        
    if(length(cases==1))
	   {
	   index<-NET[,cases]==1
	   }
       
	if(length(cases)>1)
	   {
	   index<-rowSums(NET[,cases])>=1
	   }
       
       #index<- index & !R # recovered are not really exposed
   
    E[index]<<-1
	E[cases]<<-0
	S[index]<<-0
	}

#eliminate immunity--recover sends nodes to S
recover<-function(node)
	{
	I[node]<<-0
	E[node]<<-0
	S[node]<<-1
#	R[node]<<-1
	cases<<-which(I==1)
	}



#===================================================================================
#	Simulation Setup
#	-reps: define how many time steps
#	-setup()  set S to 1, I, E, R to 0
#	-PREV: an empty vector to keep track of prevalence at each timestep
#	-par sets up a plot window. par(ask=TRUE) requires a <Enter> between plots (slows sown our movie)
#	-b : probaility infection given exposure / time step (ie rate)
#	-duration : number of days infectious
#	-Assign an index case
#	-infect it with infect()
#	-expose its neighbors with expose()
#===================================================================================

#define the length in days of the simulation
reps<-1

#set time to day 1
i<-1

#an empty vector to store the prevalence over the course of the simulation
PREV<-rep(0,reps)
PREV_A<-rep(0,reps)
PREV_B<-rep(0,reps)

#probability of infection in exposed node per time step
b <-.05

#duration of infectiousness
duration<-c(rpois(L/2,90), rpois(L/2,60))

#sets the whole populations to susceptible
setup()

#set how many index cases and draw them randomly
index<-index_case<-sample(L,5)

#infect the index cases
infect(index_case)

#expose those nodes connected to the index cases
expose()

#plot the starting conditions
par(ask=FALSE)
visualize_net()


#===================================================================================
#	Simulation loop
#	-reps: define how many time steps
#	-outputs a movie (currently commented out)
#	-plots prevalence over time
#===================================================================================

for(i in 2:reps)
	{
	random<-runif(L)
	risk<-E*b		#multiplying the risk by whether or not exposed
	infect(random<risk)
	recover(which(recoverday==i))
	expose()
    visualize_nodes()
    PREV_A[i]<-sum(cases<(L/2+1))/(L/2)
    PREV_B[i]<-sum(cases>L/2)/(L/2)
	}




plot(PREV_A,typ="l", ylab="Prevalence",xlab="Time step",ylim=c(0,1),col="blue")
lines(PREV_B, col="purple")


