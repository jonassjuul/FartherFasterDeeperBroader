package main

import (
	"fmt"
	"strconv"
	"math/rand"
    "math"
	
	"time"	
	"os"
	"log"
	"encoding/csv"
	"io"
)
func remove_element(a []int,i int) []int{
	// Remove the element at index i from a.
	a[i] = a[len(a)-1] // Copy last element to index i.
	a[len(a)-1] = 0   // Erase last element (write zero value).
	a = a[:len(a)-1]   // Truncate slice.
	return a
}

// Old functions from well-mixed case
func spread(a []string,from int) []string{
	a = append(a,strconv.Itoa(from)+"s")
	return a
}
func kill(a []string,node int) []string{
	a = append(a,strconv.Itoa(node)+"k")
	return a
}

// New functions for this particular setting...
func spread_tofrom(a []string,to int,from int) []string{
	a = append(a,strconv.Itoa(from)+"s"+strconv.Itoa(to))
	return a
}

func new_removal(num_infected int, num_SIedges int, SI_arr [][]int,cascade_string []string,network [][]int ) (int, int, [][]int,[]string) {
	var kill_node int = draw_random_integer(num_infected)
	counter := -1
	entry := 0
	// Now find out which node is actually killed...
	for counter<kill_node {

		if (len(SI_arr[entry])>0) {//&& SI_arr[entry][0]!=-1){
			if (SI_arr[entry][0]!=-1) {
				//for node_nb:=0 ; node_nb<len(SI_arr[entry]);node_nb ++{
					
				counter ++
				if (counter == kill_node) {
					// Remove the SI edges..
					num_infected --
					if (SI_arr[entry][0]!=-3) {
						num_SIedges -= len(SI_arr[entry])
					}
					SI_arr[entry] = nil
					SI_arr[entry] = append(SI_arr[entry],-1)
					cascade_string = kill(cascade_string,entry)
					break
					
				}
				//}
			}
		}
		entry ++
		// Emergency break if infectious node cannot be found...
		//if (entry < 0 || entry >= len(SI_arr)) {
		//	counter = -1
		//	kill_node --
		//	entry = 0
		//	if (kill_node < 0) {
		//		num_infected -- 
		//		fmt.Println("HERE",num_infected)
		//		break
		//
		//	}
		//}		
	}

	return num_infected, num_SIedges, SI_arr,cascade_string
}

func new_removal_IC(num_infected int, num_SIedges int, SI_arr [][]int,cascade_string []string,network [][]int ) (int, int, [][]int,[]string) {
	//var kill_node int = draw_random_integer(num_infected)

	var kill_edge int = draw_random_integer(num_SIedges)
	counter := 0
	for node := 0 ; node < len(network); node ++ {
		if (counter + len(SI_arr[node])>kill_edge ) {
			entry := kill_edge-counter 
			SI_arr[node] = remove_element(SI_arr[node],entry)
			num_SIedges -=1
			break
		} 
		counter += len(SI_arr[node])

	}


	return num_infected, num_SIedges, SI_arr,cascade_string
}


func new_infection(num_infected int, num_SIedges int, SI_arr [][]int,cascade_string []string,network [][]int ) (int, int, [][]int,[]string){
	// First, check if the infections is a seed
	var to_node int = -2
	var from_node int = -2

	if (num_infected == 0) {
		from_node = -1
		to_node = draw_random_integer(len(network))
	} else {
		var to_integer int = draw_random_integer(num_SIedges)
		to_node,from_node = find_SIedge(to_integer,SI_arr)
	}
	num_infected ++

	num_SIedges,SI_arr = update_SI_arr(to_node,num_SIedges,SI_arr,network)
	cascade_string = spread_tofrom(cascade_string,to_node, from_node)
	to_node --
	from_node --

	return num_infected,num_SIedges,SI_arr,cascade_string

}

func new_infection_IC(num_infected int, num_SIedges int, SI_arr [][]int,cascade_string []string,network [][]int ) (int, int, [][]int,[]string){
	// First, check if the infections is a seed
	var to_node int = -2
	var from_node int = -2

	if (num_infected == 0) {
		from_node = -1
		to_node = draw_random_integer(len(network))
	} else {
		var to_integer int = draw_random_integer(num_SIedges)
		to_node,from_node = find_SIedge(to_integer,SI_arr)
	}
	num_infected ++

	num_SIedges,SI_arr = update_SI_arr(to_node,num_SIedges,SI_arr,network)
	cascade_string = spread_tofrom(cascade_string,to_node, from_node)

	if (to_node != -1) {
		num_infected, num_SIedges, SI_arr,cascade_string = new_removal_IC(num_infected, num_SIedges, SI_arr,cascade_string, network)
	}

	to_node --
	from_node --

	return num_infected,num_SIedges,SI_arr,cascade_string

}



func find_SIedge (edge_num int, SI_arr [][]int) (int,int) {
	counter := -1
	to_node := -2
	from_node := -2
	for (counter < edge_num) {
		for node_from :=0 ; node_from < len(SI_arr) ; node_from ++ {
			
			if (len(SI_arr[node_from])>0 && SI_arr[node_from][0]!=-1 && SI_arr[node_from][0]!=-3) {
				for entry_to :=0 ; entry_to < len(SI_arr[node_from]) ; entry_to ++ {
					counter ++
					if (counter == edge_num) {
						from_node = node_from + 0
						to_node = SI_arr[node_from][entry_to]
						break
					}
				}
			}
			if (counter == edge_num) {
				break
			}
		}
	}
	return to_node,from_node
}

func update_SI_arr(to_node int,num_SIedges int,SI_arr [][]int,network [][]int) (int,[][]int){
	for nb_num:=0;nb_num<len(network[to_node]);nb_num ++ {
		// add new possibilities of infection
		if (len(SI_arr[network[to_node][nb_num]]) == 0) {

			// In this case, neighbour is susceptible (would have to_node as nb if infectious)
			SI_arr[to_node] = append(SI_arr[to_node],network[to_node][nb_num])
			num_SIedges ++
		} else if (SI_arr[network[to_node][nb_num]][0] != -1 && SI_arr[network[to_node][nb_num]][0] != -3)  {
			// remove out-dated possibilities of infection
			// In this case, neighbour is infectious, so must remove possbility of nb infecting to_node
			for nb_nb_num:=0;nb_nb_num<len(SI_arr[network[to_node][nb_num]]);nb_nb_num ++ {

				if (SI_arr[network[to_node][nb_num]][nb_nb_num] == to_node) {
					SI_arr[network[to_node][nb_num]] = remove_element(SI_arr[network[to_node][nb_num]],nb_nb_num)
					
					// Check if this made an infectious node with no SIedges. In that case, insert -3
					if (len(SI_arr[network[to_node][nb_num]]) == 0) {
						SI_arr[network[to_node][nb_num]] = append(SI_arr[network[to_node][nb_num]],-3) // if still infectious but no edges to spread across..
					}
					num_SIedges --
					break
				}
			}
		}
	}
	if (len(SI_arr[to_node]) == 0) {
		SI_arr[to_node] = append(SI_arr[to_node],-3) // if still infectious but no edges to spread across..


	}
	return num_SIedges,SI_arr
}


func draw_random_integer(maximum int) int {
	s1 := rand.NewSource(time.Now().UnixNano())
	r1 := rand.New(s1)
	return r1.Intn(maximum)
}
func draw_random_float() float32 {
	s1 := rand.NewSource(time.Now().UnixNano())
	r1 := rand.New(s1)
	return r1.Float32()
}
func mean_degree(network [][]int) float64 {
	number_of_edges :=0
	for node :=0 ; node<len(network) ; node ++ {
		number_of_edges+=len(network[node])
	}
	return float64(number_of_edges)/float64((len(network)))
}

func append_to_file (filename string, text []string) {
    // If the file doesn't exist, create it, or append to the file
    f, err := os.OpenFile(filename, os.O_APPEND|os.O_CREATE|os.O_WRONLY, 0644)
    if err != nil {
        log.Fatal(err)
	}
	sep := " "
	for _, line := range text {	
		if _, err := f.WriteString(line+sep); err != nil {
			log.Fatal(err)
		}

	}
	if _, err := f.WriteString("\n"); err != nil {
		log.Fatal(err)
	}	
	if err := f.Close(); err != nil {
		log.Fatal(err)
	}	
}

func divisible(number int, divisor int) bool{
	return number%divisor==0
}

func make_network_from_csv(filename string) [][]int{
	var node_number int 
	var file_content [2][]int 
	csvfile, err := os.Open(filename)
	if err != nil {
		log.Fatalln("Couldn't open the csv file", err)
	}

	// Parse the file
	r := csv.NewReader(csvfile)

	// Iterate through the records
	for {
		// Read each record from csv
		record, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatal(err)
		}
		record2,err := strconv.Atoi(record[0])
		record3,err := strconv.Atoi(record[1])

		file_content[0] = append(file_content[0],record2)
		file_content[1] = append(file_content[1],record3)


		record0,err:=strconv.ParseFloat(record[0],64)
		record1,err:=strconv.ParseFloat(record[1],64)

		if (math.Max(record0,record1)>float64(node_number)) {
			node_number = int(math.Max(record0,record1))

		}
	}
	
	// Now open again and make the array
	var network_arr = make([][]int,node_number+1)
	for i:=0; i<len(file_content[0]);i++ {
		if (file_content[0][i]!=file_content[1][i]) {
			network_arr[file_content[0][i]] = append(network_arr[file_content[0][i]],file_content[1][i])
			network_arr[file_content[1][i]] = append(network_arr[file_content[1][i]],file_content[0][i])
		}
	}
	return network_arr

}



func main() {

	// Create network
	// -------------------

	// Open the file and find the number of nodes (assuming integer numbering):

	empirical_networks :=[]string{"athletes_edges"}//{"MSU24","Cornell5","Texas84"}

	for network_id:= 0; network_id<3; network_id++ {
		


		var network_name = empirical_networks[network_id]//"MSU24"//"athletes_edges"

		// Important note: List of edges cannot list same edge twice (e.g. both 0,1 and 1,0 should not be listed.)
		
		var network [][]int = make_network_from_csv(fmt.Sprintf("%s", network_name)+".csv")//make_network_from_csv("../../../../Epidemic_descendants/ER/Edit_script_Nspreads/"+fmt.Sprintf("%s", network_name)+".csv")

		

		// Now define parameters for IC model
		const R0 float64 = 0.40

		var mean_degree float64 = mean_degree(network)
		
		var p float64 = R0/mean_degree


		// Define parameters for simulations
		const Nexp int = 10000
		var cascades_counted int = 0
		const Minimum_cascade_size int = 1
		var filename string = "Cascades_on_network_ICMODEL"+fmt.Sprintf("%s", network_name)+"_Rspread"+fmt.Sprintf("%f", R0)+"_MinimumCascadeSize"+strconv.Itoa(Minimum_cascade_size)+".txt"



		for exp:= 0 ; exp<Nexp ; exp ++ {
			fmt.Println("Doing experiment",exp,"out of",Nexp)
			var SI_arr [][]int = make([][]int,len(network))

			// first thing in each realization, choose seed
			var num_infected int 
			var num_SIedges int //= len(network[seed])

			var spread_total int = 1

			var cascade_string []string

			num_infected, num_SIedges, SI_arr, cascade_string = new_infection_IC(num_infected,num_SIedges,SI_arr,cascade_string,network)


			// Now do the following until there are no more infectious nodes.
			for (num_SIedges>0) {
				if (draw_random_float()<float32(p)) {
					num_infected, num_SIedges, SI_arr, cascade_string = new_infection_IC(num_infected,num_SIedges,SI_arr,cascade_string,network)
					spread_total ++

				} else {
					num_infected, num_SIedges, SI_arr, cascade_string = new_removal_IC(num_infected,num_SIedges,SI_arr,cascade_string,network)

				}

			}

			if (spread_total >= Minimum_cascade_size) {
				append_to_file(filename,cascade_string)
				cascades_counted += 1
				if (divisible(cascades_counted,100)){
				fmt.Println("Cascades counted:",cascades_counted)
				fmt.Printf("Last cascade had %d cases",spread_total)
				}
			}



		}
	}

}
