#!/usr/bin/awk -f
#
{
  if($2=="Number") {
    getline;
    getline;
    nat=$1;
  }
  if($5=="energy") {
    getline;
    getline;
    print $1," 1.0" >> "energies";
  }
  if($4=="gradient") {
    getline;
    for(i=1;i<=3*nat;i++) {
      getline;
      print $1," 1.0" >> "forces";
    }
  }
  if($9=="Bohr") {
    getline;
    for(i=1;i<=nat;i++) {
      getline;
      gsub(" 1 ", " H ");
      gsub(" 6 ", " C ");
      gsub(" 7 ", " N ");
      gsub(" 8 ", " O ");
      gsub(" 16 ", " S ");
      printf "%s %14.8f%14.8f%14.8f\n", $1,$2*0.529177249,$3*0.529177249,$4*0.529177249 >> "coordinates";
    }
  }
  fflush();
}
