package main

import (
	"/Helper/Import.go"
	"math/rand"
	"time"
	"fmt"
)

type mRNA struct {
	seq  string
	copy []Transcript
}

type Transcript struct {
	ribo []Ribosome
}

type Ribosome struct {
	name int

	mRNA string
	copy int
	pos int

	aSite string
	pSite string

	translating bool
}

func populateRibosomes(Number int) (map[int]Ribosome) {
	Ribosomes := make(map[int]Ribosome)

	for iN := 0; iN < Number; iN++ {
		Ribosomes[iN] = Ribosome{0, "", 0, 0, "", "", false}
	}

	return Ribosomes
}



func main() {
	rand.Seed(time.Now().UTC().UnixNano())

	// Generate maps containing a number of ribosomes as well as transcripts
	// comensurate to mRNA abundance in the cell.

	Length              := 1
	RiboMap             := populateRibosomes(Length) 
	mRNA, Counts, Names := populateTranscripts("Reference/S288C-ORFs.fa", "Reference/S288C-ORF-Counts.dat")

	fmt.Println(RiboMap)
	fmt.Println(mRNA["YEL040W"])
	fmt.Println(Counts)
	fmt.Println(Names)
	// SeqMapKeys := make([]string, 0, len(SeqMap))
	// Position   := make([]int, 0, Length)

 //    for Key := range SeqMap {
 //        SeqMapKeys = append(SeqMapKeys, Key)
 //    }
	
	// for iN := 0; iN < 10000; iN++ {
	// 	Position = make([]int, 0, Length)

	// 	for iR := 0; iR < Length; iR ++ {
	// 		RiboMap[iR], SeqMap = Simulate(RiboMap[iR], SeqMapKeys, SeqMap)
	// 		Position            = append(Position, RiboMap[iR].pos)
	// 	}

	// 	fmt.Println(RiboMap)
	// }
}