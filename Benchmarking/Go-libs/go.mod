module myapp

go 1.23.2

require github.com/consensys/gnark-crypto v0.19.2 // Fixes some security issues, including an instance of CWE-303

require (
	github.com/bits-and-blooms/bitset v1.20.0 // indirect
	github.com/consensys/bavard v0.2.1 // indirect
	github.com/dterei/gotsc v0.0.0-20160722215413-e78f872945c6
	github.com/kilic/bls12-381 v0.1.0
	github.com/mmcloughlin/addchain v0.4.0 // indirect
	golang.org/x/sys v0.30.0 // indirect
	rsc.io/tmplfunc v0.0.3 // indirect
)
