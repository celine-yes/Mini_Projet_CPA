package algorithms;

import java.awt.Point;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeSet;

public class DefaultTeam {
	
	public final static double BUDGET = 1664;
	
	
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	public Tree2D calculSteiner(ArrayList<Point> points, int edgeThreshold, ArrayList<Point> hitPoints) {
	  
		//pour definir le poids de l’arete entre deux sommets
		double[][] dist = calculShortestDist(points, edgeThreshold);
	  
		//K :Creation du graphe pondere complet 
		ArrayList<Edge> graphComplet = pointsToGraph(points, dist, hitPoints);
		 
		//T0 :arbre couvrant à partir de kruskal sur K
		ArrayList<Edge> kruskal = Kruskal(graphComplet, hitPoints);
		  
		int [][] path = calculShortestPaths(points, edgeThreshold);
		//Construction de graphe H en remplacant toute arete uv par un plus court chemin entre u et v
		ArrayList<Edge> H = getShortestPath(kruskal, path, dist, points);
		  
		//reconstruction de l'arbre couvrant la plus petite 
		ArrayList<Edge> kruskal2 = Kruskal(H, points);
		return edgesToTree(kruskal2, kruskal2.get(0).pointA);
	}
	
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
    public Tree2D calculSteinerBudget(ArrayList<Point> points, int edgeThreshold, ArrayList<Point> hitPoints) {
    	
    	double[][] dist = calculShortestDist(points, edgeThreshold); //pour determiner les poids des aretes dans le graphe
		ArrayList<Edge> K = pointsToGraph(points, dist, hitPoints); //graphe pondéré complet
        ArrayList<Edge> kruskal = Kruskal(K, points); //'algorithme de Kruskal sur le graphe complet pour trouver un arbre couvrant minimum
        
        ArrayList<Edge> greedy = new ArrayList<>();
        
        //pour facilite la recherche de l'arête avec la plus petite distance pour chaque point
        Map<Point, TreeSet<Edge>> pointTreeSet = pointToTreeSetEdges(kruskal);
        
        //une liste pour suivre les points déjà ajoutés à l'arbre
        ArrayList<Point> addedPoints = new ArrayList<>();
        addedPoints.add(hitPoints.get(0)); // maison-mère
        
        double totalCost = 0;
        while (!pointTreeSet.isEmpty() && totalCost <= BUDGET) {
        	
        	//trouver l'arête avec la distance minimale connectée à un point déjà ajouté
            Edge minDistEdge = null;
            double minDist = Double.POSITIVE_INFINITY;
	        
	        for (Point p : addedPoints) {
	        	//ensemble trié des arêtes associées à ce point
	            TreeSet<Edge> edges = pointTreeSet.get(p);
	            if (edges != null && !edges.isEmpty()) {
	            	//l'arete ayant la plus petite distance
	                Edge edge = edges.first();
	                //on cherche le minimum parmi toutes les aretes min
	                if (edge.dist < minDist) {
	                    minDist = edge.dist;
	                    minDistEdge = edge;
	                }
	            }
	        }
            
            Point pointA = minDistEdge.pointA;
			Point pointB = minDistEdge.pointB;

			// Ajouter les deux extremites (points) de cette arete aux points ajoutes pour la recherche de min plus tard
			if (!addedPoints.contains(pointA))
				addedPoints.add(pointA);
			pointTreeSet.get(pointA).remove(minDistEdge);
			
			if (!addedPoints.contains(pointB))
				addedPoints.add(pointB);
			pointTreeSet.get(pointB).remove(minDistEdge);
			
			if (totalCost + minDistEdge.dist <= BUDGET) {
	            greedy.add(minDistEdge);
	            totalCost += minDistEdge.dist;
			}else{
				break;
			}
        }

        int[][] chemin = calculShortestPaths(points, edgeThreshold);
        //remplacer toute arete uv par un plus court chemin entre u et v
        ArrayList<Edge> gs = getShortestPath(greedy, chemin, dist, points);
        return edgesToTree(gs, hitPoints.get(0));
    }
    
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
	  //algo floyd-warshall pour calculer les chemins les plus courts entre tous les points
    public int[][] calculShortestPaths(ArrayList<Point> points, int edgeThreshold) {
		  
		int n = points.size();
	    int[][] paths=new int[n][n];
	    double[][] dist= new double[n][n];
	    
	    for (int i=0;i<paths.length;i++) {
	    	for (int j=0;j<paths.length;j++) {
	    		paths[i][j]=i; 
	    	}
	    }
	    
	    // Initialisation des distances avec les distances directes entre les points
	    for (int i = 0; i < n; i++) {
	        for (int j = 0; j < n; j++) {
	            Point p1 = points.get(i);
	            Point p2 = points.get(j);
	            if (i==j) {
	            	// Distance nulle pour le meme point
	            	dist[i][i]=0; 
	            	continue;
	            }
	            
	            double distance = p1.distance(p2);
	            if (distance <= edgeThreshold) {
	                dist[i][j] = distance;
	            } else {
	                dist[i][j] = Double.POSITIVE_INFINITY;
	            }
	            paths[i][j]=j;
	        }
	    }
	    
	    for (int k=0; k<n ; k++) {
	    	for(int i=0; i<n ; i++) {
	    		for(int j=0; j<n ; j++) {
	    			if (dist[i][k] + dist[k][j] < dist[i][j]) {
	                    dist[i][j] = dist[i][k] + dist[k][j];
	                    // l’indice du prochain sommet dans un plus court chemin de points i a j
	                    paths[i][j]=paths[i][k];
	                }
	    		}
	    	}
	    }
	    return paths;
	  }
	 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
	  
	//algo floyd-warshall pour calculer les distances de chemin le plus court entre tous les points
    public double[][] calculShortestDist(ArrayList<Point> points, int edgeThreshold) {
		  
	    double[][] dist= new double[points.size()][points.size()];
	    
	    // Initialisation des distances avec les distances directes entre les points
	    for (int i = 0; i < points.size(); i++) {
	        for (int j = 0; j < points.size(); j++) {
	            Point p1 = points.get(i);
	            Point p2 = points.get(j);
	            if (i==j) {
	            	// Distance nulle pour le meme point
	            	dist[i][i]=0; 
	            	continue;
	            }
	            
	            double distance = p1.distance(p2);
	            if (distance <= edgeThreshold) {
	                dist[i][j] = distance;
	            } else {
	                dist[i][j] = Double.POSITIVE_INFINITY;
	            }
	        }
	    }
	    
	    for (int k=0; k<points.size() ; k++) {
	    	for(int i=0; i<points.size() ; i++) {
	    		for(int j=0; j<points.size() ; j++) {
	    			if (dist[i][k] + dist[k][j] < dist[i][j]) {
	                    dist[i][j] = dist[i][k] + dist[k][j];
	                }
	    		}
	    	}
	    }
	    return dist;
	  }
	 
	  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
	  
	  public ArrayList<Edge> pointsToGraph (ArrayList<Point> points, double[][] dist, ArrayList<Point> hitPoints){ 
		  
		  ArrayList<Edge> edges = new ArrayList<Edge>();
		  
		  // Création de toutes les combinaison d'arêtes entre les points de hitPoints = graphe pondere complet
		  for (int i = 0; i < hitPoints.size()-1; i++) {
		      for (int j = i+1 ; j < hitPoints.size(); j++) {
		    	  //les index pour recuperer apres le chemin le plus court
		          int indexI = points.indexOf(hitPoints.get(i));
		          int indexJ = points.indexOf(hitPoints.get(j));
		          edges.add(new Edge(points.get(indexI), 
					        		 points.get(indexJ), 
					        		 indexI, 
					        		 indexJ,
		        		  			dist[indexI][indexJ])); //poids = la longeur du plus court chemin entre point i et j
		      }
		  }
		  return edges;
	  }
	  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
	  //renvoie toutes les points du plus court chemin entre 2 points 
	  public ArrayList<Edge> getShortestPath (ArrayList<Edge> edges, int[][] path, double[][] dist,  ArrayList<Point> points){
		  ArrayList<Edge> shortestPath = new ArrayList<Edge>();
		  
		  for(int i=0; i<edges.size(); i++) {
			  //deux extremites d'arete
			  int indexStart = edges.get(i).indexA; 
			  int end = edges.get(i).indexB;
			  //l'indice suivant sur le chemin
			  int indNextPoint = path[indexStart][end];
			  
			  while(indNextPoint != indexStart) { 
				  shortestPath.add(new Edge(points.get(indexStart), 
						  					points.get(indNextPoint), 
						  					indexStart, 
						  					indNextPoint,
						  					dist[indexStart][indNextPoint]));
				  indexStart = indNextPoint;
				  //parce que paths[i][j] renvoie l’indice k du prochain sommet dans un plus court chemin
				  indNextPoint = path[indexStart][end];
				  
			  }
		  }
		  return shortestPath;

	  }
	  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	  
	  
	  //renvoie des aretes pour construire l'arbre couvrant avec la plus petite longuer total des aretes
	  public ArrayList<Edge> Kruskal(ArrayList<Edge> graph, ArrayList<Point> points) {
		  
		  //triées par ordre croissant de poids pour les explorer en premier 
		  Collections.sort(graph);

		  ArrayList<Edge> kruskal = new ArrayList<Edge>();
		  Edge current=null;
		  LabelPoint forest = new LabelPoint(points);
		  
		  while (graph.size()!=0) {
		    current = graph.remove(0);
		    //pour eviter des cycles
		    if (forest.getLabel(current.pointA)!=forest.getLabel(current.pointB)) {
		      kruskal.add(current);
		      //renommage des etiquettes des points pour suivre les composantes connexes du graphe
		      forest.setLabel(forest.getLabel(current.pointA),forest.getLabel(current.pointB));
		    }
		  }
		  return kruskal;
		}
	  
	  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	  private Map<Point, TreeSet<Edge>> pointToTreeSetEdges(ArrayList<Edge> edges) {
	    	
			//associer chaque point à l'ensemble des arêtes connectées à ce point 
			//ces aretes sont trié par leur poids
			//cette structure est cree pour optimiser le temps de calcul
			Map<Point, TreeSet<Edge>> edgesMap = new HashMap<>();
			
			for (Edge edge : edges) {
			    Point pointA = edge.pointA;
			    Point pointB = edge.pointB;
			    if (!edgesMap.containsKey(pointA)) {
			    	edgesMap.put(pointA, new TreeSet<>());
			    }
			    //le trie se fait automatiquement dans TreeSet avec la fonction compareTo de la classe Edge
			    edgesMap.get(pointA).add(edge);
			    if (!edgesMap.containsKey(pointB)) {
			    	edgesMap.put(pointB, new TreeSet<>());
			    }
			    edgesMap.get(pointB).add(edge);
			}
			return edgesMap;
	    }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
	  //renvoie un arbre construit récursivement sur les sous-arbres de la racine spécifiée
	  @SuppressWarnings("unchecked")
	private Tree2D edgesToTree(ArrayList<Edge> edges, Point racine) {
		  	
		    //les arêtes restantes après avoir traité les arêtes qui ont racine comme point de départ ou d'arrivée
			ArrayList<Edge> reste = new ArrayList<Edge>();
			// sous-arbre de la racine
			ArrayList<Point> racinesSousArbres = new ArrayList<Point>();
			Edge actuelle = null;
			
			while (edges.size() != 0) {
				actuelle = edges.remove(0);
					//si point depart de l'arete = racine
				if (actuelle.pointA.equals(racine)) {
					racinesSousArbres.add(actuelle.pointB);
				} else {
					//si point arrive = racine
					if (actuelle.pointB.equals(racine)) {
						racinesSousArbres.add(actuelle.pointA);
					} else {
						reste.add(actuelle);
					}
				}
			}
			//appels récursifs pour compléter l'arbre principal
			ArrayList<Tree2D> sousArbres = new ArrayList<Tree2D>();
			for (Point racineSousArbre : racinesSousArbres) {
				sousArbres.add(edgesToTree((ArrayList<Edge>) reste.clone(), racineSousArbre));
			}
			
			return new Tree2D(racine, sousArbres);
		}
	  


}
