package main

import (
	"math/rand"
	"time"
	"fmt"

	// Required for populateTranscripts
	"strings"
	"strconv"
	"bufio"
	"os"

	// Required for simulate
	"sort"
)

type mRNA struct {
	seq     string
	copy    []Transcript
	count   int
	riboTot int
}

type Transcript struct {
	ribo  []Ribosome
}

type Ribosome struct {
	name int

	mRNA string
	copy int
	pos int

	aSite string
	pSite string

	trans bool
}

func populateRibosomes(Number int) (map[int]Ribosome) {
	Ribosomes := make(map[int]Ribosome)

	for iN := 0; iN < Number; iN++ {
		Ribosomes[iN] = Ribosome{iN, "", 0, -1, "", "", false}
	}

	return Ribosomes
}

func buildTAIMap(TAIFile string) (map[string]float64) {

	// Scan a list of TAI values and if 

	InFile, Err := os.Open(TAIFile)
	if Err != nil { panic(Err) }
	defer InFile.Close()

	Scanner := bufio.NewScanner(InFile)
	TAIMap  := make(map[string]float64)

	for Scanner.Scan() {
		Split            := strings.Split(Scanner.Text(), " ")
		TAI, _           := strconv.ParseFloat(Split[1], 64)
		TAIMap[Split[0]]  = TAI
	}

	return TAIMap
}

func populateTranscripts(FASTA string, Counts string) (map[string]mRNA, []float64, []string) {

	// Scan through a file of transcript counts for genes across an organism. ThiStep helps
	// the simulation better reflect the actual biology as often times a few hundred genes
	// help characterize the transcriptional landscape of the organism.
	//
	// A map of gene name to transcript counts iStep generated: Count

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
	// Scan through a file of transcript counts for genes across an organism. ThiStep helps
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
					mRNAs[Name] = mRNA{Sequence, make([]Transcript, Count[Name]), Count[Name], 0}
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

func simStep(R Ribosome, RAll map[int]Ribosome, C []float64, N []string, M map[string]mRNA, T map[string]float64) (Ribosome, map[string]mRNA) {
	if (R.trans == false) {
		if (rand.Float64() > 0.975) {

			// Randomly initiate on one of the genes. The spacings of the random number floats
			// between 0 to 1 reflects the abundance of the different transcripts.
		
			Check  := rand.Float64()
			Idx    := sort.Search(len(C), func(i int) bool { return C[i] >= Check }) - 1

			PmRNA  := N[Idx]
			PCopy  := rand.Intn(M[PmRNA].count)
			PTrans := M[PmRNA].copy[PCopy]
			Seq    := M[PmRNA].seq[0:6]

			// Once we've defined the (P)rospective mRNA and Copy to initiate on, we check to see
			// that there isn't already a ribosome (the most recently initiated one), in the way.
			// We define that to be a clearance of at least 18 nucleotides. If thiStep condition is
			// met, then the ribosome can initiate with no issue.

			if (len(PTrans.ribo) == 0) {
				R = Ribosome{R.name, PmRNA, PCopy, 0, Seq[3:6], Seq[0:3], true}
				M[R.mRNA].copy[R.copy].ribo = append(M[R.mRNA].copy[R.copy].ribo, R)

				var tMap      = M[R.mRNA]
				tMap.riboTot += 1
				M[R.mRNA]     = tMap
			} else {
				if (RAll[PTrans.ribo[len(PTrans.ribo) - 1].name].pos > 18) {
					R = Ribosome{R.name, PmRNA, PCopy, 0, Seq[3:6], Seq[0:3], true}
					M[R.mRNA].copy[R.copy].ribo = append(M[R.mRNA].copy[R.copy].ribo, R)
					
					var tMap      = M[R.mRNA]
					tMap.riboTot += 1
					M[R.mRNA]     = tMap
				}
			}			
		}
	} else {
		if (R.pos + 7 > len(M[R.mRNA].seq)) {
			if (rand.Float64() > 0.9) {
				M[R.mRNA].copy[R.copy].ribo = M[R.mRNA].copy[R.copy].ribo[1:]
				R      = Ribosome{R.name, "", 0, -1, "", "", false}
				
				var tMap      = M[R.mRNA]
				tMap.riboTot -= 1
				M[R.mRNA]     = tMap
			}
		} else {
			if (rand.Float64() < T[R.aSite]) { 
				RList := M[R.mRNA].copy[R.copy].ribo
				Pos   := 0

				for iX := 0; iX < len(RList); iX++ {
					if (R.name == RList[iX].name) {
						Pos = iX
						break
					}
				}

				if (Pos != 0) {
					Latter := RAll[RList[Pos - 1].name].pos
					Former := RAll[RList[Pos].name].pos

					if (Latter - Former) > 18 {
						R.pos   = R.pos + 3
						R.pSite = M[R.mRNA].seq[R.pos:(R.pos + 3)]
	 					R.aSite = M[R.mRNA].seq[(R.pos + 3):(R.pos + 6)]

	 					M[R.mRNA].copy[R.copy].ribo[Pos] = R
					}
				} else {
					R.pos   = R.pos + 3
					R.pSite = M[R.mRNA].seq[R.pos:(R.pos + 3)]
	 				R.aSite = M[R.mRNA].seq[(R.pos + 3):(R.pos + 6)]

 					M[R.mRNA].copy[R.copy].ribo[Pos] = R
				}		
			}
		}
	}
	return R, M
}

func main() {
	rand.Seed(time.Now().UTC().UnixNano())

	// Generate the data structures for all the ribosomes as well as all the transcripts in the "cell".
	// 
	// The mRNA stucture contains the sequence of the mRNA as well as a map to all individual transcripts.
	// These transcripts in turn will have a slice of ribosomes which can be used to detect "collisions."
	// The final information stored in the mRNA structure iStep the total counts of mRNA to make rand.Intn
	// calls nice and easy.
	//
	// The slice "Counts" contains floats ranging from 0 to 1 with non-even spacing. The size of the intervals
	// represents the relative amount of transcripts in the transcript pool so that initiation iStep commensurate
	// to transcript abundance.
	//
	// Finally, the slice "Name" contains the strings corresponding to the names of the transcripts above.

	Start := time.Now()

	Length              := 100000
	RiboMap             := populateRibosomes(Length) 
	TAIMap              := buildTAIMap("Reference/TAI.dat")
	mRNA, Counts, Names := populateTranscripts("Reference/S288C-ORFs.fa", "Reference/S288C-ORF-Counts.dat")

	mRNACounts, E := os.Create("./Reads.dat")
	if E != nil { fmt.Println("Check permissions? Can't make data file!"); os.Exit(1) }
	defer mRNACounts.Close()

	mRNACountsWriter    := bufio.NewWriter(mRNACounts)

	mRNACountsWriter.WriteString("Step" + "\t")

	for iName := range(Names) {
		mRNACountsWriter.WriteString(Names[iName] + "\t")
	}

	mRNACountsWriter.WriteString("\n")

	for iStep := 1; iStep < 25001; iStep++ {

		for iR := 0; iR < Length; iR++ {
			RiboMap[iR], mRNA = simStep(RiboMap[iR], RiboMap, Counts, Names, mRNA, TAIMap)
		}

		if (iStep % 250 == 0) {
			mRNACountsWriter.WriteString(strconv.Itoa(iStep)+ "\t")

    		for iName := range(Names) {
    			mRNACountsWriter.WriteString(strconv.Itoa(mRNA[Names[iName]].riboTot) + "\t")
    		}

    		mRNACountsWriter.WriteString("\n")
		}
	}

	mRNACountsWriter.Flush()
	fmt.Println(time.Since(Start))
}
