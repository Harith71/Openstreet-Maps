// application.cpp <Starter Code>
// Harith Patel
//
// University of Illinois at Chicago
// CS 251: Spring 2022
// Project #7 - Openstreet Maps
//

// CREATIVE COMPONENT: finds a buulding that is closer to a person that is
// more tired and prints paths to the building for each person.
//
// TO RUN: enter c when asked to run creative component or application
//         enter building for person 1 and person 2
//         enter 1 or 2 when asked which person is more tired

// References:
// TinyXML: https://github.com/leethomason/tinyxml2
// OpenStreetMap: https://www.openstreetmap.org
// OpenStreetMap docs:
//   https://wiki.openstreetmap.org/wiki/Main_Page
//   https://wiki.openstreetmap.org/wiki/Map_Features
//   https://wiki.openstreetmap.org/wiki/Node
//   https://wiki.openstreetmap.org/wiki/Way
//   https://wiki.openstreetmap.org/wiki/Relation
//

#include <iostream>
#include <iomanip> /*setprecision*/
#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <queue>
#include <limits>

#include "tinyxml2.h"
#include "dist.h"
#include "osm.h"
#include "graph.h"

using namespace std;
using namespace tinyxml2;

const double INF = numeric_limits<double>::max();

// Compares edge weights for priority queue in Dijkstra's algorithm
class prioritize {
  public:
  bool operator()(pair<long long, double> &p1, 
                  pair<long long, double> &p2) const {
    return p1.second > p2.second;
  }
};

// PrintPath:
//
// prints the shortest path from starting vertex to the ending vertex
void printPath(map<long long, long long> pred, long long endV) {
  vector<long long> path;
  long long currV = endV;
  // loop until starting verext is reached
  while (currV != 0) {
    path.push_back(currV);
    currV = pred[currV];
  }
  // iterate vector backwards since the path is added backwards into array.
  cout << "Path: ";
  for (int i = path.size() - 1; i >= 0; i--) {
    if (i == 0) {
      cout << path[i];
      continue;
    }
    cout << path[i] << "->";
  }
  cout << endl;
}

// printBuildingInfo:
//
// prints building names and coordinates of each building
void printBuildingInfo(BuildingInfo b1, BuildingInfo b2, BuildingInfo dest) {
  cout << endl
       << "Person 1's point:" << endl;
  cout << " " << b1.Fullname << endl;
  cout << " (" << b1.Coords.Lat << ", " << b1.Coords.Lon << ")" << endl;
  cout << "Person 2's point:" << endl;
  cout << " " << b2.Fullname << endl;
  cout << " (" << b2.Coords.Lat << ", " << b2.Coords.Lon << ")" << endl;
  cout << "Destination Building:" << endl;
  cout << " " << dest.Fullname << endl;
  cout << " (" << dest.Coords.Lat << ", " << dest.Coords.Lon << ")" << endl;
}

// printNewDest:
//
// prints info of the new closest building found
void printNewDest(BuildingInfo b, long long n,
                  map<long long, Coordinates> &Nodes) {
  cout << endl
       << "New destination building:" << endl;
  cout << " " << b.Fullname << endl;
  cout << " (" << b.Coords.Lat << ", " << b.Coords.Lon << ")" << endl;
  cout << "Nearest destination node:" << endl;
  cout << " " << n << endl;
  cout << " (" << Nodes[n].Lat << ", " << Nodes[n].Lon << ")" << endl;
}

// printNodeInfo:
//
// prints nodeID and coordinates of each node
void printNodeInfo(long long n1, long long n2, long long dN,
                   map<long long, Coordinates> &Nodes) {
  cout << endl
       << "Nearest P1 node:" << endl;
  cout << " " << n1 << endl;
  cout << " (" << Nodes[n1].Lat << ", " << Nodes[n1].Lon << ")" << endl;
  cout << "Nearest P2 node:" << endl;
  cout << " " << n2 << endl;
  cout << " (" << Nodes[n2].Lat << ", " << Nodes[n2].Lon << ")" << endl;
  cout << "Nearest destination node:" << endl;
  cout << " " << dN << endl;
  cout << " (" << Nodes[dN].Lat << ", " << Nodes[dN].Lon << ")" << endl;
}

// DijkstraShortestPath:
//
// finds the shortest distance from staring vertex to every other
// reachable vertex in the graph and updtaes predecessors and distances
// map with the path to each vertex.
void DijkstraShortestPath(long long s, graph<long long, double> &G,
                          map<long long, long long> &predecessors,
                          map<long long, double> &distances) {
  priority_queue<pair<long long, double>,
                 vector<pair<long long, double>>, prioritize>
      unvisitedQueue;
  vector<long long> verts = G.getVertices();

  // set default weight and predecessor values for each vertex in the maps
  for (long long v : verts) {
    predecessors[v] = 0;
    distances[v] = INF;
  }

  // change distance of starting vertex and add to the front of priority queue
  distances[s] = 0;
  unvisitedQueue.push(make_pair(s, 0));

  while (!unvisitedQueue.empty()) {
    long long currV = unvisitedQueue.top().first;
    double currW = unvisitedQueue.top().second;
    unvisitedQueue.pop();

    // check if current vertex is not unreachable.
    if (distances[currV] == INF)
      break;

    set<long long> adjVert = G.neighbors(currV);
    for (long long adjV : adjVert) {
      // find distance to vertex.
      double edgeW;
      if (!G.getWeight(currV, adjV, edgeW))
        continue;
      double altPathW = currW + edgeW;

      // if new path is shorter than path in map update it to the new path
      if (altPathW < distances[adjV]) {
        distances[adjV] = altPathW;
        predecessors[adjV] = currV;
        // add the adjacent vertex and its weight to priority queue.
        unvisitedQueue.push(make_pair(adjV, altPathW));
      }
    }
  }
}

// searchBuilding:
//
// searches the vector with buildings to see if given name string matches
// any building if there is a match the buulding struct is returned.
// Otherwise defualt building struct is returned.
BuildingInfo searchBuilding(string name, vector<BuildingInfo> &Buildings) {
  BuildingInfo b;

  for (int i = 0; i < Buildings.size(); i++) {
    if (Buildings[i].Abbrev == name)
      return Buildings[i];

    // check building name contians given search string
    if (Buildings[i].Fullname.find(name) != string::npos)
      return Buildings[i];
  }
  return b;
}

// nearestBuild:
//
// finds the nearest building to given coordinates if building is in
// the unreachable buildings vector it can't be the closes buulding.
// The closes building to the coordinates given is returned.
BuildingInfo nearestBuild(Coordinates c, vector<BuildingInfo> &Buildings,
                          vector<string> &unrBuild) {
  double min = INF;
  int pos = -1;
  double dist;
  for (int i = 0; i < Buildings.size(); i++) {
    int found = 0;
    BuildingInfo b = Buildings[i];
    // check if building is not in unreachable building vector
    for (int i = 0; i < unrBuild.size(); i++) {
      if (unrBuild[i] == b.Fullname) {
        found = 1;
        break;
      }
    }

    if (found == 1)
      continue;
    dist = distBetween2Points(c.Lat, c.Lon, b.Coords.Lat, b.Coords.Lon);

    if (dist < min) {
      min = dist;
      pos = i;
    }
  }
  return Buildings[pos];
}


// nearestNode:
//
// finds the nearest node to a building that is given. The nodeID
// of the nearest node is returned.
long long nearestNode(BuildingInfo b, map<long long, Coordinates> &Nodes,
                      vector<FootwayInfo> &Footways) {
  double min = INF;
  double dist;
  long long node;
  // loop through every node in each footway
  for (auto fw : Footways) {
    for (int i = 0; i < fw.Nodes.size(); i++) {
      // get coordinates of node in footway
      long long n1 = fw.Nodes[i];
      Coordinates c1 = Nodes[n1];
      dist = distBetween2Points(c1.Lat, c1.Lon, b.Coords.Lat, b.Coords.Lon);

      if (dist < min) {
        min = dist;
        node = n1;
      }
    }
  }
  return node;
}

// application:
//
// finds a destination between 2 people and prints the shortest path
// to the destination building from the location of each person.
void application(
    map<long long, Coordinates> &Nodes, vector<FootwayInfo> &Footways,
    vector<BuildingInfo> &Buildings, graph<long long, double> &G) {
  string person1Building, person2Building;
  Coordinates midP;
  BuildingInfo dest, b1, b2;
  long long b1N, b2N, destN;

  cout << endl;
  cout << "Enter person 1's building (partial name or abbreviation), or #> ";
  getline(cin, person1Building);

  while (person1Building != "#") {
    map<long long, long long> predecessors1, predecessors2;
    map<long long, double> distances1, distances2;
    vector<string> unrBuild;
    int attempts = 0;
    bool found = false;
    cout << "Enter person 2's building (partial name or abbreviation)> ";
    getline(cin, person2Building);

    b1 = searchBuilding(person1Building, Buildings);
    b2 = searchBuilding(person2Building, Buildings);

    while (!found) {
      if (b1.Fullname.size() == 0) {
        cout << "Person 1's building not found" << endl;
        break;
      }

      if (b2.Fullname.size() == 0) {
        cout << "Person 2's building not found" << endl;
        break;
      }

      midP = centerBetween2Points(b1.Coords.Lat, b1.Coords.Lon,
                                  b2.Coords.Lat, b2.Coords.Lon);
      dest = nearestBuild(midP, Buildings, unrBuild);
      b1N = nearestNode(b1, Nodes, Footways);
      b2N = nearestNode(b2, Nodes, Footways);
      destN = nearestNode(dest, Nodes, Footways);

      // check if it is the first try trying to
      // find a destination building and path
      if (attempts == 0) {
        printBuildingInfo(b1, b2, dest);
        printNodeInfo(b1N, b2N, destN, Nodes);
      } else {
        printNewDest(dest, destN, Nodes);
      }
      DijkstraShortestPath(b1N, G, predecessors1, distances1);
      DijkstraShortestPath(b2N, G, predecessors2, distances2);

      // check if there is a path between person 1 and person 2
      if (distances1[b2N] >= INF) {
        cout << "Sorry, destination unreachable." << endl;
        break;
      }

      // check if each person has a path to the destination building
      if (distances1[destN] >= INF || distances2[destN] >= INF) {
        cout << endl
             << "At least one person was unable to reach the destination building.";
        cout << " Finding next closest building..." << endl;
        attempts = 1;
        unrBuild.push_back(dest.Fullname);
        continue;
      }

      cout << endl
           << "Person 1's distance to dest: " << distances1[destN] << " miles" << endl;
      printPath(predecessors1, destN);
      cout << endl
           << "Person 2's distance to dest: " << distances2[destN] << " miles" << endl;
      printPath(predecessors2, destN);
      break;
    }

    cout << endl;
    cout << "Enter person 1's building (partial name or abbreviation), or #> ";
    getline(cin, person1Building);
  }
}

// creative:
//
// finds a destination building that is closer to the more tired person and
// prints the shortest path to the destination building from each person.
void creative(
    map<long long, Coordinates> &Nodes, vector<FootwayInfo> &Footways,
    vector<BuildingInfo> &Buildings, graph<long long, double> &G) {
    string person1Building, person2Building;
    Coordinates midP;
    BuildingInfo dest, b1, b2;
    long long b1N, b2N, destN;

    cout << endl;
    cout << "Enter person 1's building (partial name or abbreviation), or #> ";
    getline(cin, person1Building);

    while (person1Building != "#") {
      map<long long, long long> predecessors1, predecessors2;
      map<long long, double> distances1, distances2;
      vector<string> unrBuild;
      int attempts = 0;
      bool found = false;
      string person;
      cout << "Enter person 2's building (partial name or abbreviation)> ";
      getline(cin, person2Building);

      b1 = searchBuilding(person1Building, Buildings);
      b2 = searchBuilding(person2Building, Buildings);

      while (!found) {
        if (b1.Fullname.size() == 0) {
          cout << "Person 1's building not found" << endl;
          break;
        }

        if (b2.Fullname.size() == 0) {
          cout << "Person 2's building not found" << endl;
          break;
        }
        cout << "Enter the person number that is more tired(1 or 2): ";
        getline(cin, person);
        midP = centerBetween2Points(b1.Coords.Lat, b1.Coords.Lon,
                                    b2.Coords.Lat, b2.Coords.Lon);
        cout << midP.Lat << " " << midP.Lon << endl;
        // find a position that is closer to the person that is more tired
        if (person == "1") {
          midP = centerBetween2Points(midP.Lat, midP.Lon, b1.Coords.Lat, b1.Coords.Lon);
          midP = centerBetween2Points(midP.Lat, midP.Lon, b1.Coords.Lat, b1.Coords.Lon);
        } else {
          midP = centerBetween2Points(midP.Lat, midP.Lon, b2.Coords.Lat, b2.Coords.Lon);
          midP = centerBetween2Points(midP.Lat, midP.Lon, b2.Coords.Lat, b2.Coords.Lon);
        }
        
        cout << midP.Lat << " " << midP.Lon << endl;
        dest = nearestBuild(midP, Buildings, unrBuild);
        b1N = nearestNode(b1, Nodes, Footways);
        b2N = nearestNode(b2, Nodes, Footways);
        destN = nearestNode(dest, Nodes, Footways);

        // check if it is the first try trying to
        // find a destination building and path
        if (attempts == 0) {
          printBuildingInfo(b1, b2, dest);
          printNodeInfo(b1N, b2N, destN, Nodes);
        } else {
          printNewDest(dest, destN, Nodes);
        }
        DijkstraShortestPath(b1N, G, predecessors1, distances1);
        DijkstraShortestPath(b2N, G, predecessors2, distances2);

        if (distances1[b2N] >= INF) {
          cout << "Sorry, destination unreachable." << endl;
          break;
        }

        if (distances1[destN] >= INF || distances2[destN] >= INF) {
          cout << endl
              << "At least one person was unable to reach the destination building.";
          cout << " Finding next closest building..." << endl;
          attempts = 1;
          unrBuild.push_back(dest.Fullname);
          continue;
        }

        cout << endl
            << "Person 1's distance to dest: " << distances1[destN] << " miles" << endl;
        printPath(predecessors1, destN);
        cout << endl
            << "Person 2's distance to dest: " << distances2[destN] << " miles" << endl;
        printPath(predecessors2, destN);
        break;
      }

      cout << endl;
      cout << "Enter person 1's building (partial name or abbreviation), or #> ";
      getline(cin, person1Building);
    }
}

int main() {
  // maps a Node ID to it's coordinates (lat, lon)
  map<long long, Coordinates> Nodes;
  // info about each footway, in no particular order
  vector<FootwayInfo> Footways;
  // info about each building, in no particular order
  vector<BuildingInfo> Buildings;
  XMLDocument xmldoc;

  cout << "** Navigating UIC open street map **" << endl;
  cout << endl;
  cout << std::setprecision(8);

  string def_filename = "map.osm";
  string filename;

  cout << "Enter map filename> ";
  getline(cin, filename);

  if (filename == "") {
    filename = def_filename;
  }

  //
  // Load XML-based map file
  //
  if (!LoadOpenStreetMap(filename, xmldoc)) {
    cout << "**Error: unable to load open street map." << endl;
    cout << endl;
    return 0;
  }

  //
  // Read the nodes, which are the various known positions on the map:
  //
  int nodeCount = ReadMapNodes(xmldoc, Nodes);

  //
  // Read the footways, which are the walking paths:
  //
  int footwayCount = ReadFootways(xmldoc, Footways);

  //
  // Read the university buildings:
  //
  int buildingCount = ReadUniversityBuildings(xmldoc, Nodes, Buildings);

  //
  // Stats
  //
  assert(nodeCount == (int)Nodes.size());
  assert(footwayCount == (int)Footways.size());
  assert(buildingCount == (int)Buildings.size());

  cout << endl;
  cout << "# of nodes: " << Nodes.size() << endl;
  cout << "# of footways: " << Footways.size() << endl;
  cout << "# of buildings: " << Buildings.size() << endl;

  graph<long long, double> G;

  // add nodes into graps as vertex
  for (auto n : Nodes)
    G.addVertex(n.first);

  // add egdes into graph
  for (auto fw : Footways) {
    for (int i = 0; i < fw.Nodes.size() - 1; i++) {
      long long n1 = fw.Nodes[i];
      long long n2 = fw.Nodes[i + 1];
      Coordinates c1 = Nodes[n1];
      Coordinates c2 = Nodes[n2];
      // get distance between node and the node after it
      double distance = distBetween2Points(c1.Lat, c1.Lon, c2.Lat, c2.Lon);
      // add edge from one to other and vice versa so it isn't a directed graph
      G.addEdge(n1, n2, distance);
      G.addEdge(n2, n1, distance);
    }
  }

  cout << "# of vertices: " << G.NumVertices() << endl;
  cout << "# of edges: " << G.NumEdges() << endl;
  cout << endl;

  //
  // Menu
  //
  string userInput;
  cout << "Enter \"a\" for the standard application or "
       << "\"c\" for the creative component application> ";
  getline(cin, userInput);
  if (userInput == "a") {
    application(Nodes, Footways, Buildings, G);
  } else if (userInput == "c") {
    creative(Nodes, Footways, Buildings, G);
  }
  //
  // done:
  //
  cout << "** Done **" << endl;
  return 0;
}
