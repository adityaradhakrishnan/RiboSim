package main

import (
	"bufio"
	"fmt"
	"os"
)

type Ribosome struct {
	name, aSite, pSite string
	pos, posRefStart, posRefStop int
	translating bool
}

type mRNA struct {
	name, sequence string
}

func populateTranscripts(FASTA string, Counts string) {
	InFile, Err  := os.Open(FASTA)
	if Err != nil { panic(Err) }
	defer InFile.Close()

	Scanner := bufio.NewScanner(InFile)

	fmt.Println(Counts)

	for Scanner.Scan() {
		Line := Scanner.Text()

		if (Line[:1] == ">") {
			fmt.Println(Line)
		}
	}

}

func main() {
	R := Ribosome{"1", "GAA", "ATG", 0, -15, 12, true}

	populateTranscripts("MG1655v3-ORFs.fa", "CAGATC-DhOE.fastq")

	// Print default zeroed value ribosome
	fmt.Println("Default ribosome is: ", R)
}