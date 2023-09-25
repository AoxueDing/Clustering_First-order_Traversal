package cmsc420_s23;

import java.sql.Wrapper;
import java.util.ArrayList;
import java.util.HashMap;

public class FarthestFirst<LPoint extends LabeledPoint2D> {
	ClusterAssignment cluster;
	WtLeftHeap heap;
	HashMap<String, WtLeftHeap.Locator> map;
	ArrayList<LPoint> traversalList;
	public FarthestFirst(int rebuildOffset, Rectangle2D bbox, LPoint startCenter) throws Exception {
		cluster = new ClusterAssignment(rebuildOffset, bbox, startCenter);
		heap = new WtLeftHeap();
		map = new HashMap<String, WtLeftHeap.Locator>();
		traversalList = new ArrayList<>();
		traversalList.add(startCenter);
	}
	public void addSite(LPoint site) throws Exception {
		cluster.addSite(site);
		AssignedPair pair = new AssignedPair();
		pair.site = site;
		pair.center = cluster.closestCenter(site);
		pair.distanceSq = site.getPoint2D().distanceSq(pair.center.getPoint2D());
		WtLeftHeap.Locator loc = heap.insert(pair.distanceSq, pair);
		map.put(site.getLabel(), loc);
	}
	public LPoint extractNext() throws Exception{
		if(cluster.sitesSize() == 0){
			return null;
		}
		AssignedPair pair = (AssignedPair) heap.extract();
		Point2D site =  pair.site.getPoint2D();
		cluster.kdTree.delete(site);
		ArrayList<LPoint> sites = cluster.addCenterWithReturn(pair.site);
		traversalList.add((LPoint) pair.site);
		for(LPoint curr: sites){
			WtLeftHeap.Locator loc = map.get(curr.getLabel());
			double newDis = curr.getPoint2D().distanceSq(site);
			heap.updateKey(loc, newDis);
		}
		return (LPoint) pair.site;
	}
	public int sitesSize() { /* ... */ return cluster.sitesSize(); }
	public int traversalSize() { /* ... */ return traversalList.size(); }
	public void clear()throws Exception {
		cluster.clear();
		map.clear();
		traversalList.clear();
	}
	public ArrayList<String> listKdWithCenters() { /* ... */ return cluster.listKdWithCenters(); }
	public ArrayList<LPoint> getTraversal() { /* ... */ return traversalList; }
	public ArrayList<String> listCenters() { /* ... */ return cluster.listCenters(); }
	public ArrayList<String> listAssignments() { /* ... */ return cluster.listAssignments(); }
}
