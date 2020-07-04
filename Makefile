#
# Makefile for Microbial Community Analysis
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# written by Conrad Shyu (conradshyu@hotmail.com)
#
# Author's comments:
# ------------------
# Institute of Bioinformatics and Computational Biology
# University of Idaho, Moscow, ID 83844
#
# revised on September 3, 2008
# revised on March 12, 2014
#
all: erpa ispar pat pspa trflp

erpa:
	g++ -O3 -I. bitvector.cpp cmdparam.cpp seqdb.cpp erpa.cpp -o erpa -lpthread
ispar:
	g++ -O3 -I. bitvector.cpp cmdparam.cpp seqdb.cpp ispar.cpp -o ispar -lpthread
pat:
	g++ -O3 -I. bitvector.cpp cmdparam.cpp seqdb.cpp pat.cpp -o pat -lpthread
pspa:
	g++ -O3 -I. bitvector.cpp cmdparam.cpp seqdb.cpp pspa.cpp -o pspa -lpthread
trflp:
	g++ -O3 -I. bitvector.cpp cmdparam.cpp seqdb.cpp trflp.cpp -o trflp -lpthread

clean:
	rm -f erpa ispar pat pspa trflp
