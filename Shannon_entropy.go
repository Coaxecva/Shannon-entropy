package main

import (
	"fmt"
	"os"
	"bufio"
	"bytes"
	"io/ioutil"
	"math"
	"strconv"
)

func Shannon_entropy(seq []byte, winSize int, step int) float64 {
	m := make(map[string]int)
	str := string(seq)
	num := ((len(seq)-winSize)/step)+1
	for i := 0; i<num*step; i+=step {
		_, ok := m[str[i:i+winSize]]
		if !ok {
			m[str[i:i+winSize]] = 1
		} else {
			m[str[i:i+winSize]] += 1
		}
	}
	
	s := 0.0
	for _, value := range m {
		p := (float64(value) / float64(len(seq) - winSize + 1))
		s += p*math.Log2(p)
	}
	return (-s)
}

func main() {
	if len(os.Args) != 3 {
      panic("must provide sequence file.")
   	}
   
   	seq := ReadSequence(os.Args[1])
   	k, _ := strconv.Atoi(os.Args[2])
   	fmt.Println(Shannon_entropy(seq, k, 1))
}

func ReadSequence(file string) []byte {
   f, err := os.Open(file)
   if err != nil {
      panic(err)
   }
   defer f.Close()
   byte_array := make([]byte, 0)
   Ns := []byte("N")
   None := []byte("")
   if file[len(file)-6:] == ".fasta" || file[len(file)-3:] == ".fa" {
      scanner := bufio.NewScanner(f)
      for scanner.Scan() {
         line := scanner.Bytes()
         if len(line)>0 && line[0] != '>' {
            byte_array = append(byte_array, bytes.Replace(bytes.Trim(line,"\n\r "), Ns, None, -1)...)
         }
      }
   } else {
      byte_array, err = ioutil.ReadFile(file)
      if err != nil {
         panic(err)
      }
   }
   return byte_array
}