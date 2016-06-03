package main

import (
	"strings"
	"strconv"
	"bufio"
	"os"
)

func populateTranscripts(FASTA string, Counts string) (map[string]mRNA, []float64, []string) {

	// Scan through a file of transcript counts for genes across an organism. This helps
	// the simulation better reflect the actual biology as often times a few hundred genes
	// help characterize the transcriptional landscape of the organism.
	//
	// A map of gene name to transcript counts is generated: Count

	InFile, Err  := os.Open(Counts)
	if Err != nil { panic(Err) }
	defer InFile.Close()

	Scanner := bufio.NewScanner(InFile)
	Count   := make(map[string]int)

	for Scanner.Scan() {
		Split            := strings.Split(Scanner.Text(), "\t")
		Length, _        := strconv.Atoi(Split[1])
		Count[Split[0]]   = Length
	}

	// Scan through the FASTA file and populate a map of mRNAs
	// Scan through a file of transcript counts for genes across an organism. This helps
	// the simulation better reflect the actual biology as often times a few hundred genes
	// help characterize the transcriptional landscape of the organism.

	InFile, Err  = os.Open(FASTA)
	if Err != nil { panic(Err) }
	defer InFile.Close()

	Scanner   = bufio.NewScanner(InFile)
	mRNAs    := make(map[string]mRNA)
	mRNAName := []string{}
	mRNACnt  := []float64{0}
	Sequence := ""
	Name     := ""
	Counter  := 0.0

	for Scanner.Scan() {
		Line := Scanner.Text()

		if (Line[:1] == ">") {		
			if (len(Sequence) > 0) {		
				if (Count[Name] > 0) {	
					Counter    += float64(Count[Name])
					mRNAs[Name] = mRNA{Sequence, make([]Transcript, Count[Name])}
					mRNAName    = append(mRNAName, Name)
					mRNACnt     = append(mRNACnt, mRNACnt[len(mRNACnt) - 1] + float64(Count[Name]))
				}
			}
			Name          = Line[1:] 
			Sequence      = ""
		} else {
			Sequence = Sequence + Line[:len(Line)]
		}
	}

	for iN := 0; iN < len(mRNAName); iN++ {
		mRNACnt[iN] = mRNACnt[iN]/Counter
	}

	return mRNAs, append(mRNACnt[:len(mRNACnt)-1], 1), mRNAName
}