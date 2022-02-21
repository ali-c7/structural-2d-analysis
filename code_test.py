import numpy as np
import math

class Member:
    '''
    Fields:
        Length (Nat),
        Theta (Nat),
        Global DOF (listof Nat),
        Local DOF (listof Nat),
        DOF Restraints (listof Str)
    '''

    def __init__(self, number1):
        '''
        Constructor: Create a Member object by
        calling Member(start_x, end_x, start_y, end_y, theta, global_dof, local_dof, restraints)

        Effects: Mutates Self

        __init__: Member Nat Nat Nat Nat Nat (listof Nat) (listof Nat) (listof Str) -> None
        '''

        self.label = st.markdown("Member " + str(number1))
        self.name = "Member #" + str(number1)

        self.col1, self.col2, self.col3, self.col4, self.col5, self.col6, self.col7, self.col8 = st.columns(8)
        with self.col1:
            self.start_node = float(st.text_input('Start', value=0, key="node_start_"+str(number1)))

        with self.col2:
            self.end_node = float(st.text_input('End', value=0, key="node_end_" + str(number1)))

        with self.col3:
            self.x_start_dof = float(st.text_input("x DOF start", value=0, key="dof_start_x" + str(number1)))

        with self.col4:
            self.y_start_dof = float(st.text_input("y DOF  start", value=0, key="dof_start_y"+ str(number1)))

        with self.col5:
            self.z_start_dof = float(st.text_input("z DOF start", value=0, key="dof_start_z"+ str(number1)))

        with self.col6:
            self.x_end_dof = float(st.text_input("x DOF end", value=0, key="dof_end_x"+ str(number1)))

        with self.col7:
            self.y_end_dof = float(st.text_input("y DOF end", value=0, key="dof_end_y"+ str(number1)))

        with self.col8:
            self.z_end_dof = float(st.text_input("z DOF end ", value=0, key="dof_end_z"+ str(number1)))

        self.start_x = nodes[int(self.start_node) - 1].x
        self.start_y = nodes[int(self.start_node) - 1].y
        self.end_x = nodes[int(self.end_node) - 1].x
        self.end_y = nodes[int(self.end_node) - 1].y
        self.id = number1
        self.member_length = math.sqrt((self.end_x-self.start_x) ** 2 + (self.end_y-self.start_y) ** 2)
        self.start_node_dof = [self.x_start_dof, self.y_start_dof, self.z_start_dof]
        self.end_node_dof = [self.x_end_dof, self.y_end_dof, self.z_end_dof]

        if self.start_x == self.end_x:
                if self.start_y < self.end_y:
                    self.theta = 90  #* np.pi / 180
                else:
                    self.theta = -90
        else:
            self.theta = math.atan((self.end_y-self.start_y)/(self.end_x-self.start_x)) * 180 / np.pi
        self.label = st.sidebar.markdown('m' + str(number1) + ":" + f"  L = {str(self.member_length)}m," + f" Angle = {str(self.theta)}" )

    def __repr__(self):
        s = "Member #{0.id}"
        return s.format(self)

    def member_matrix(self):
        AE = 6300000    # AE
        EI = 105000  # EI
        if self.start_x == self.end_x:
                if self.start_y < self.end_y:
                    self.theta = 90
                else:
                    self.theta = -90
        else:
            self.theta = math.atan((self.end_y-self.start_y)/(self.end_x-self.start_x)) * 180 / np.pi
        member_stiffness_matrix = np.array([
            [((AE / self.member_length) * (np.cos(self.theta * np.pi / 180)) ** 2) + (
                    (12 * EI) / (self.member_length ** 3)) * (np.sin(self.theta * np.pi / 180)) ** 2,
             ((AE / self.member_length) - (12 * EI / (self.member_length ** 3))) * (
                     np.cos(self.theta * np.pi / 180) * np.sin(self.theta * np.pi / 180)),
             (-6 * EI / (self.member_length ** 2)) * np.sin(self.theta * np.pi / 180), -(
                    ((AE / self.member_length) * (np.cos(self.theta * np.pi / 180)) ** 2) + (
                    (12 * EI) / (self.member_length ** 3)) * (np.sin(self.theta * np.pi / 180)) ** 2),
             -((AE / self.member_length) - (12 * EI / (self.member_length ** 3))) * (
                     np.cos(self.theta * np.pi / 180) * np.sin(self.theta * np.pi / 180)),
             (-6 * EI / (self.member_length ** 2)) * np.sin(self.theta * np.pi / 180)],
            [((AE /self.member_length) - (12 * EI / (self.member_length ** 3))) * (
                    np.cos(self.theta * np.pi / 180) * np.sin(self.theta * np.pi / 180)),
             ((AE / self.member_length) * (np.sin(self.theta * np.pi / 180)) ** 2) + (12 * EI / (self.member_length ** 3)) *
             (np.cos(self.theta * np.pi / 180)) ** 2,
             (6 * EI / (self.member_length ** 2)) * np.cos(self.theta * np.pi / 180), -(
                    ((AE / self.member_length) - (12 * EI / (self.member_length ** 3))) * (
                    np.cos(self.theta * np.pi / 180) * np.sin(self.theta * np.pi / 180))), -1 * (
                     ((AE / self.member_length) * (np.sin(self.theta * np.pi / 180)) ** 2) + (
                     (12 * EI) / (self.member_length ** 3)) * (np.cos(self.theta * np.pi / 180)) ** 2),
             (6 * EI / (self.member_length ** 2)) * np.cos(self.theta * np.pi / 180)],
            [(-6 * EI / (self.member_length ** 2)) * np.sin(self.theta * np.pi / 180),
             (6 * EI / (self.member_length ** 2)) * np.cos(self.theta * np.pi / 180), (4 * EI / self.member_length),
             (6 * EI / (self.member_length ** 2)) * np.sin(self.theta * np.pi / 180),
             (-6 * EI / (self.member_length ** 2)) * np.cos(self.theta * np.pi / 180), (2 * EI / self.member_length)],
            [-(((AE / self.member_length) * (np.cos(self.theta * np.pi / 180)) ** 2) + (
                    (12 * EI) / (self.member_length ** 3)) * (np.sin(self.theta * np.pi / 180)) ** 2), -(
                    ((AE / self.member_length) - (12 * EI / (self.member_length ** 3))) * (
                    np.cos(self.theta * np.pi / 180) * np.sin(self.theta * np.pi / 180))),
             (6 * EI / (self.member_length ** 2)) * np.sin(self.theta * np.pi / 180), (
                     ((AE / self.member_length) * (np.cos(self.theta * np.pi / 180)) ** 2) + (
                     (12 * EI) / (self.member_length ** 3)) * (np.sin(self.theta * np.pi / 180)) ** 2),
             ((AE / self.member_length) - (12 * EI / (self.member_length ** 3))) * (
                     np.cos(self.theta * np.pi / 180) * np.sin(self.theta * np.pi / 180)),
             (6 * EI / (self.member_length ** 2)) * np.sin(self.theta * np.pi / 180)],
            [-((AE / self.member_length) - (12 * EI / (self.member_length ** 3))) * (
                    np.cos(self.theta * np.pi / 180) * np.sin(self.theta * np.pi / 180)), -1 * (
                     ((AE / self.member_length) * (np.sin(self.theta * np.pi / 180)) ** 2) + (
                     (12 * EI) / (self.member_length ** 3)) * (np.cos(self.theta * np.pi / 180)) ** 2),
             (-6 * EI / (self.member_length ** 2)) * np.cos(self.theta * np.pi / 180),
             ((AE / self.member_length) - (12 * EI / (self.member_length ** 3))) * (
                     np.cos(self.theta * np.pi / 180) * np.sin(self.theta * np.pi / 180)),
             ((AE / self.member_length) * (np.sin(self.theta * np.pi / 180)) ** 2) + (12 * EI / (self.member_length ** 3)) * (
                 np.cos(self.theta * np.pi / 180)) ** 2,
             (-6 * EI / (self.member_length ** 2)) * np.cos(self.theta * np.pi / 180)],
            [(-6 * EI / (self.member_length ** 2)) * np.sin(self.theta * np.pi / 180),
             (6 * EI / (self.member_length ** 2)) * np.cos(self.theta * np.pi / 180), 2 * EI / self.member_length,
             (6 * EI / (self.member_length ** 2)) * np.sin(self.theta * np.pi / 180),
             (-6 * EI / (self.member_length ** 2)) * np.cos(self.theta * np.pi / 180), (4 * EI / self.member_length)]])

        return member_stiffness_matrix


def global_stiffness_matrix_fwd_pass(lom):
    # consumes list of members
    # checks if beginning/ending dof's are the same
    gsm = np.zeros(((len(lom)+1)*3, (len(lom)+1)*3))  # need to fix, global stiffness matrix changes
    end_dof = get_end_dof(lom)
    start_dof = get_start_dof(lom)
    for i in range(len(lom)):
        for j in range(len(lom)):
            print("Comparing DOF")
            print(lom[i].start_node_dof, lom[j].end_node_dof)
            if lom[i].start_node_dof == lom[j].end_node_dof and lom[i].start_node_dof in end_dof and end_dof.count(lom[i].start_node_dof) > 1 and j == end_dof.index(lom[i].start_node_dof):
                start = int(lom[i].start_node_dof[0])
                stop = int(lom[i].start_node_dof[-1])
                gsm[start-1:stop, start-1:stop] += lom[i].member_matrix()[0:3, 0:3] + lom[j].member_matrix()[3:,3:]
                print("Condition 1 ran")
            elif lom[i].start_node_dof == lom[j].end_node_dof and lom[i].start_node_dof in end_dof and end_dof.count(lom[i].start_node_dof) > 1:
                start = int(lom[i].start_node_dof[0])
                stop = int(lom[i].start_node_dof[-1])
                gsm[start-1:stop, start-1:stop] +=  lom[j].member_matrix()[3:,3:]
                print("Condition 2 ran")
            elif lom[i].start_node_dof == lom[j].end_node_dof and lom[i].start_node_dof in end_dof and end_dof.count(lom[i].start_node_dof) == 1:
                start = int(lom[i].start_node_dof[0])
                stop = int(lom[i].start_node_dof[-1])
                print("lom[i]")
                print(lom[i].member_matrix())
                print("lom[j]")
                print(lom[j].member_matrix())
                gsm[start - 1:stop, start - 1:stop] += np.add(lom[i].member_matrix()[0:3, 0:3], lom[j].member_matrix()[3:, 3:])
                print("Condition 3 ran")
            elif lom[i].start_node_dof != lom[j].end_node_dof and j == len(lom) - 1 and lom[i].start_node_dof not in end_dof:
                start = int(lom[i].start_node_dof[0])
                stop = int(lom[i].start_node_dof[-1])
                gsm[start - 1:stop, start - 1:stop] += lom[i].member_matrix()[0:3, 0:3]
                print("Condition 4 ran")
            print(gsm)
            print("----------------------------")
    start = int(max(end_dof)[0] - 1)
    stop = int(max(end_dof)[-1])
    gsm[start:stop, start:stop] += lom[end_dof.index(max(end_dof))].member_matrix()[3:, 3:] # have to fix, it's hardcoded
    return gsm