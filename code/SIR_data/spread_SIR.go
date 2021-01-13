package main

import (
	"fmt"
	"strconv"
    "math/rand"
	"time"	
	"os"
	"log"
)
func remove_element(a []int,i int) []int{
	// Remove the element at index i from a.
	a[i] = a[len(a)-1] // Copy last element to index i.
	a[len(a)-1] = 0   // Erase last element (write zero value).
	a = a[:len(a)-1]   // Truncate slice.
	return a
}
func spread(a []string,from int) []string{
	a = append(a,strconv.Itoa(from)+"s")
	return a
}
func kill(a []string,node int) []string{
	a = append(a,strconv.Itoa(node)+"k")
	return a
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




func main() {



	// SIR parameters..
	const R_death float32 = 1
	const R_spread float32 = 0.85

	const p_spread float32 = R_spread/(R_death+R_spread)

	// Cascade parameters..
	const N_cascades = 30000
	const Minimum_cascade_size = 50
	var cascades_counted int = 0

	// Define file name
	var filename string = "Cascades_Rspread"+fmt.Sprintf("%f", R_spread)+"_MinimumCascadeSize"+strconv.Itoa(Minimum_cascade_size)+".txt"

	for cascades_counted < N_cascades {
		// For each cascade, keep track of..
		var events []string
		var living_nodes []int
		living_nodes = append(living_nodes,0)

		var spread_total int = 1

		// dynamics..
		for len(living_nodes) > 0 {	

			var node_affected int = draw_random_integer(len(living_nodes))//living_nodes[draw_random_integer(len(living_nodes))]
			if (draw_random_float()<p_spread) {

				events = spread(events,living_nodes[node_affected])
				living_nodes = append(living_nodes,spread_total)
				spread_total += 1

			} else {
				events = kill(events,living_nodes[node_affected])
				living_nodes = remove_element(living_nodes,node_affected)

			}
		}
		if (spread_total >= Minimum_cascade_size) {
			append_to_file(filename,events)
			cascades_counted += 1
			if (divisible(cascades_counted,100)){
			fmt.Println("Cascades counted:",cascades_counted)
			fmt.Printf("Last cascade had %d cases",spread_total)
			}
		}


	}
}
