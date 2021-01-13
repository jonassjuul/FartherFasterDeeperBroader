// You can edit this code!
// Click here and start typing.
package main

import (
	"fmt"
	"strconv"
	"math/rand"
	"math"
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
func draw_random_float() float64 {
	s1 := rand.NewSource(time.Now().UnixNano())
	r1 := rand.New(s1)
	return r1.Float64()
}

func factorial (integer float64) float64{

	var res float64=  1
	var dummy_index float64 = 1
	for dummy_index <= integer {
		res *= dummy_index
		dummy_index += 1
	}
	
	return res
}

func Poisson_dist (mean_dist float64, integer float64) float64 {

	return math.Pow(mean_dist,integer)*math.Exp(-mean_dist)/factorial(integer)
}

func draw_random_Poisson(mean_dist float64) float64 {
	// Draw a random number between 0 and 1
	s1 := rand.NewSource(time.Now().UnixNano())
	r1 := rand.New(s1)
	rand_float := r1.Float64()	

	// Find first integer, i, where sum(P(x),x<=i) < random number
	var sum_var float64 = 0
	var integer float64 = -1
	for sum_var < rand_float {
		integer  += 1
		sum_var += Poisson_dist(mean_dist,integer)
	}

	return integer
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



	// IC parameters..

	const R0 float64 = 0.90

	// Cascade parameters..
	const N_cascades = 60000
	const Minimum_cascade_size = 1
	var cascades_counted int = 0

	// Define file name
	var filename string = "ICCascades_Rspread"+fmt.Sprintf("%f", R0)+"_MinimumCascadeSize"+strconv.Itoa(Minimum_cascade_size)+".txt"

	for cascades_counted < N_cascades {

		var events []string
		

		var living_node_min int = 0
		var living_node_max int = 0

		var spread_total int = 1

		// dynamics..
		for living_node_max-living_node_min >= 0 {	

			var num_of_children int = int(draw_random_Poisson(R0))

			
			for child:=1; child <= num_of_children; child ++ {
				events = spread(events,living_node_min)
				living_node_max += 1
				spread_total += 1


			}

			events = kill(events,living_node_min)
			living_node_min +=1

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
