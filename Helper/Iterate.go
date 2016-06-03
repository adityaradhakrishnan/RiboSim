
func Simulate(R Ribosome, S []string, SM map[string]mRNA) (Ribosome, map[string]mRNA) {
	if R.translating == false {
		if (rand.Float64() > 0.975) {

			// Randomly initiate on a mRNA. More abundant transcripts are
			// more likely to be initiated on. Once you've chosen an mRNA

			R.mRNA  = S[rand.Intn(len(S))]
			Copies := len(SM[R.mRNA].tCopy)
			R.copy  = rand.Intn(Copies)


			SM[R.mRNA].transCopy[R.copy] += 1

			R.pSite = SM[R.mRNA].seq[R.pos:(R.pos + 3)]
			R.aSite = SM[R.mRNA].seq[(R.pos + 3):(R.pos + 6)] 

			R.translating = true
		}

	} else {
		if (R.pos + 7 > len(SM[R.mRNA].seq)) {
			if (rand.Float64() > 0.9) {
				R = Ribosome{"", 0, 0, "", "", false}
			}
		} else {
			if (rand.Float64() > 0.75) {
				R.pos   = R.pos + 3
				R.pSite = SM[R.mRNA].seq[R.pos:(R.pos + 3)]
				R.aSite = SM[R.mRNA].seq[(R.pos + 3):(R.pos + 6)] 
			}
		}
	}

	return R, SM
}
