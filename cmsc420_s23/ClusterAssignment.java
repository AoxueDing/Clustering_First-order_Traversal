package cmsc420_s23;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;

public class ClusterAssignment<LPoint extends LabeledPoint2D> {
	LPoint startCenter;
	SMkdTree<LPoint> kdTree; // storage for the sites
	private ArrayList<LPoint> centers = new ArrayList<>(); // storage for the centers
	public ClusterAssignment(int rebuildOffset, Rectangle2D bbox, LPoint startCenter) throws Exception {
		kdTree = new SMkdTree<>(rebuildOffset, bbox, startCenter);
		kdTree.addCenter(startCenter, null);
		centers.add(startCenter);
		this.startCenter = startCenter;

	}
	public void addSite(LPoint site) throws Exception {
		kdTree.insert(site);
	}
	public void deleteSite(LPoint site) throws Exception {
		kdTree.delete( site.getPoint2D());
	}
	public void addCenter(LPoint center) throws Exception {
		//kdTree.addCenter(center,null);
		addCenterWithReturn(center);
	}

	public ArrayList<LPoint> addCenterWithReturn(LPoint center) throws Exception {
		centers.add(center);
		return kdTree.addCenter(center, null);
	}

	public int sitesSize() { /* ... */ return kdTree.size(); }
	public int centersSize() { /* ... */ return centers.size(); }
	public void clear() throws Exception {
		kdTree.clear();
		kdTree.addCenter(startCenter, null);
		for (int i = centers.size() - 1; i >= 0; i--) {
			centers.remove(i);
		}
		centers.add(startCenter);
	}
	public ArrayList<String> listKdWithCenters() {

			ArrayList<String> ans = kdTree.listWithCenters();
			return ans;
			}
	public ArrayList<String> listCenters() {
		Collections.sort(centers, (p1, p2) -> p1.getLabel().compareTo(p2.getLabel()));
		ArrayList<String> ans = new ArrayList<String>();
		ans.add(getCenters(centers));
		return ans;
	}

	private String getCenters(ArrayList<LPoint> centers) {
		StringBuilder sb = new StringBuilder();
		for (LPoint curr : centers) {
			sb.append(curr.toString()).append("\n  ");
		}
		if (sb.length() > 0) {
			sb.delete(sb.length() - 3, sb.length()); // remove the last "\n  "
		}
		return sb.toString();
	}


	public ArrayList<String> listAssignments() {
		ArrayList<AssignedPair> pairs = kdTree.listAssignments();
		ArrayList<String> ans = assignmentString(pairs);
		return ans;

	}

	private ArrayList<String> assignmentString(ArrayList<AssignedPair> pairs) {
		ArrayList<String> result = new ArrayList<>();
		for (AssignedPair pair : pairs) {
			String pairStr = "[" + pair.site.getLabel()
					+ "->" + pair.center.getLabel()+ "]"
					+ " distSq = " + pair.distanceSq;
			result.add(pairStr);
		}
		return result;
	}

	public LPoint closestCenter(LPoint site) {
		return kdTree.closestCenter(site);
	}

    public void deleteCenter(LPoint ap) {
    }

//	// based on center return the nearest sites
//	public <LPoint extends LabeledPoint2D> ArrayList<LPoint> getSites(LabeledPoint2D site, Rectangle2D cell, ArrayList<LPoint> list) {
//		return kdTree.getSites(site, null,null);
//	}
}
