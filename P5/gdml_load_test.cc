// P5 standalone gate 0: local G4 11.4.2 parses the exported sPHENIX GDML.
#include <G4GDMLParser.hh>
#include <G4LogicalVolume.hh>
#include <G4VPhysicalVolume.hh>
#include <iostream>
int main(int argc, char** argv)
{
  G4GDMLParser p;
  p.Read(argv[1], false);
  G4VPhysicalVolume* w = p.GetWorldVolume();
  std::cout << "GDMLTEST world: " << w->GetName()
            << " daughters: " << w->GetLogicalVolume()->GetNoDaughters() << std::endl;
  return 0;
}
