use crate::matrix::Matrix;
use num_traits::{One, Zero};
use std::ops::{AddAssign, Div, Mul, MulAssign, Neg};

impl<T> Matrix<T>
where
    T: Zero + One + PartialEq,
{
    fn is_in_zeilenstufenform(&self) -> bool {
        let mut blocked_columns = 0;
        for row in 0..self.rows {
            let pivot_column = self.cells[row]
                .iter()
                .take_while(|cell| cell.is_zero())
                .count();
            if pivot_column < blocked_columns {
                return false;
            } else {
                blocked_columns = pivot_column + 1;
            }
        }

        true
    }

    fn is_in_strict_zeilenstufenform(&self) -> bool {
        let mut blocked_columns = 0;
        let mut used_columns = vec![false; self.columns];
        for row in 0..self.rows {
            let pivot_column = self.cells[row]
                .iter()
                .take_while(|cell| cell.is_zero())
                .count();
            if pivot_column < blocked_columns {
                return false;
            } else {
                blocked_columns = pivot_column + 1;

                if pivot_column < self.columns {
                    if used_columns[pivot_column] {
                        return false;
                    }

                    for column in 0..self.columns {
                        let cell = &self[(row, column)];
                        if !cell.is_zero() {
                            used_columns[column] = true;
                        }
                    }
                }
            }
        }

        true
    }

    fn is_in_exact_solution_form(&self) -> bool {
        let mut one_column = 0;
        for row in 0..self.rows {
            for column in 0..self.columns {
                let cell = &self[(row, column)];
                if one_column == self.columns - 1 {
                    if !cell.is_zero() {
                        return false;
                    } else {
                        // pass
                    }
                } else if column == self.columns - 1 {
                    // pass
                } else if column == one_column {
                    if !cell.is_one() {
                        return false;
                    } else {
                        // pass
                    }
                } else {
                    if !cell.is_zero() {
                        return false;
                    } else {
                        // pass
                    }
                }
            }
            one_column += 1;
        }

        true
    }
}

impl<T> Matrix<T>
where
    T: Default
        + Zero
        + One
        + PartialEq
        + Mul<Output = T>
        + MulAssign
        + AddAssign
        + Div<Output = T>
        + Neg<Output = T>
        + Clone,
{
    fn multiply_row(&mut self, row: usize, factor: T) {
        for column in 0..self.columns {
            self[(row, column)] *= factor.clone();
        }
    }

    fn add_multiplied_row(&mut self, origin_row: usize, factor: T, target_row: usize) {
        for column in 0..self.columns {
            let origin_row_cell = self[(origin_row, column)].clone();
            self[(target_row, column)] += origin_row_cell * factor.clone();
        }
    }

    fn turn_into_zeilenstufenform(&mut self) {
        let mut finished_rows = 0;
        loop {
            if finished_rows == self.rows || finished_rows + 1 == self.rows {
                assert!(self.is_in_zeilenstufenform());
                return;
            }

            let mut pivot_row_column = None;
            'column_loop: for column in 0..self.columns {
                for row in finished_rows..self.rows {
                    let cell = &self[(row, column)];
                    if !cell.is_zero() {
                        pivot_row_column = Some((row, column));
                        break 'column_loop;
                    }
                }
            }

            let Some((pivot_row, pivot_column)) = pivot_row_column else {
                assert!(self.is_in_zeilenstufenform());
                return;
            };

            self.swap_rows(finished_rows, pivot_row);
            let pivot_row = finished_rows;

            let pivot_cell = self[(pivot_row, pivot_column)].clone();
            for target_row in (pivot_row + 1)..self.rows {
                let cell_below_pivot = self[(target_row, pivot_column)].clone();
                let factor = -(cell_below_pivot / pivot_cell.clone());
                self.add_multiplied_row(pivot_row, factor, target_row);
                assert!(self[(target_row, pivot_column)].is_zero());
            }

            finished_rows += 1;
        }
    }

    fn turn_into_strict_zeilenstufenform(&mut self) {
        self.turn_into_zeilenstufenform();

        let mut used_columns = vec![false; self.columns];
        for row in 0..self.rows {
            let pivot_column = self.cells[row]
                .iter()
                .take_while(|cell| cell.is_zero())
                .count();
            if pivot_column < self.columns {
                if used_columns[pivot_column] {
                    let pivot_cell = self[(row, pivot_column)].clone();
                    for target_row in 0..row {
                        let target_row_cell = self[(target_row, pivot_column)].clone();
                        let factor = -(target_row_cell / pivot_cell.clone());
                        self.add_multiplied_row(row, factor, target_row);
                        assert!(self[(target_row, pivot_column)].is_zero());
                    }
                    used_columns[pivot_column] = false;
                }

                for column in 0..self.columns {
                    let cell = &self[(row, column)];
                    if !cell.is_zero() {
                        used_columns[column] = true;
                    }
                }
            }
        }

        assert!(self.is_in_strict_zeilenstufenform());
    }

    pub fn solve_as_system_of_linear_equations(&mut self) -> LinearEquationSystemSolutions<T> {
        self.turn_into_strict_zeilenstufenform();

        'next_row: for row in 0..self.rows {
            for column in 0..self.columns - 1 {
                if !self[(row, column)].is_zero() {
                    continue 'next_row;
                }
            }
            if !self[(row, self.columns - 1)].is_zero() {
                return LinearEquationSystemSolutions::NoSolution;
            }
        }

        let mut offset = Matrix::new(self.columns - 1, 1);
        let mut parameter_matrix = Matrix::new(self.columns - 1, 0);

        for row in 0..self.rows {
            let involved_variables = (0..self.columns - 1)
                .filter(|column| !self[(row, *column)].is_zero())
                .collect::<Vec<_>>();

            if involved_variables.len() == 0 {
                break;
            } else {
                let first_variable_column = involved_variables[0];
                let first_variable_cell = self[(row, first_variable_column)].clone();
                let rhs_cell = self[(row, self.columns - 1)].clone();
                offset[(first_variable_column, 0)] = rhs_cell / first_variable_cell;

                if involved_variables.len() == 1 {
                    continue;
                }

                /*

                involved variables: 0, 1, 2
                               row: 3, 4, 5
                                    3x + 4y + 5z = 0

                a := 3x
                b := 4y

                x = 1/3 * a
                y = 1/4 * b
                z = -3/5 * x - 4/5 * y = -1/5 * a - 1/5 * b

                (  1/3   0  )   ( a )
                (   0   1/4 ) * ( b )
                ( -1/5 -1/5 )

                 */

                let (&last_variable_column, other_variables) =
                    involved_variables.split_last().unwrap();

                let first_parameter_column = parameter_matrix.columns;
                parameter_matrix.add_columns(other_variables.len());

                let last_variable_cell = self[(row, last_variable_column)].clone();

                for (parameter_offset, &variable_column) in other_variables.iter().enumerate() {
                    let variable_cell = self[(row, variable_column)].clone();
                    let parameter_column = first_parameter_column + parameter_offset;
                    parameter_matrix[(variable_column, parameter_column)] =
                        T::one() / variable_cell;
                    parameter_matrix[(last_variable_column, parameter_column)] =
                        -(T::one() / last_variable_cell.clone());
                }
            }
        }

        if parameter_matrix.columns == 0 {
            LinearEquationSystemSolutions::ExactSolution(offset)
        } else {
            LinearEquationSystemSolutions::ManySolutions {
                offset,
                parameter_matrix,
            }
        }
    }
}

#[derive(Debug)]
pub enum LinearEquationSystemSolutions<T> {
    NoSolution,
    ExactSolution(Matrix<T>),
    ManySolutions {
        offset: Matrix<T>,
        parameter_matrix: Matrix<T>,
    },
}

#[cfg(test)]
mod tests {
    use crate::matrix::Matrix;

    #[test]
    fn test_multiply() {
        let a = Matrix::from_cells([[1, 2], [3, 4]]);
        let b = Matrix::from_cells([[5, 6], [7, 8]]);
        let c = Matrix::from_cells([[9, 10], [11, 12]]);
        println!("{:?}", (a.clone() * b.clone()) * c.clone());
        println!("{:?}", a * (b * c));
    }

    #[test]
    fn test_1() {
        let matrix = Matrix::from_cells([[0, 1, 2], [1, 0, 0], [0, 0, 0]]);
        assert_eq!(matrix.is_in_zeilenstufenform(), false);
        assert_eq!(matrix.is_in_exact_solution_form(), false);
    }

    #[test]
    fn test_2() {
        let matrix = Matrix::from_cells([[0, 1, 2], [0, 1, 1], [0, 0, 0]]);
        assert_eq!(matrix.is_in_zeilenstufenform(), false);
        assert_eq!(matrix.is_in_exact_solution_form(), false);
    }

    #[test]
    fn test_3() {
        let matrix = Matrix::from_cells([[1, 2, -1], [0, 0, -1], [0, 0, 0]]);
        assert_eq!(matrix.is_in_zeilenstufenform(), true);
        assert_eq!(matrix.is_in_strict_zeilenstufenform(), false);
        assert_eq!(matrix.is_in_exact_solution_form(), false);
    }

    #[test]
    fn test_4() {
        let matrix = Matrix::from_cells([[1, 2, 0], [0, 0, -1], [0, 0, 0]]);
        assert_eq!(matrix.is_in_strict_zeilenstufenform(), true);
        assert_eq!(matrix.is_in_exact_solution_form(), false);
    }

    #[test]
    fn test_5() {
        let matrix = Matrix::from_cells([[1, 0, 3], [0, 1, -1], [0, 0, 0]]);
        assert_eq!(matrix.is_in_strict_zeilenstufenform(), true);
        assert_eq!(matrix.is_in_exact_solution_form(), true);
    }

    #[test]
    fn test_6() {
        let matrix = Matrix::from_cells([[-3.0, 2.0, 0.0], [3.0, -2.0, 0.0], [1.0, 1.0, 1.0]]);
        assert_eq!(matrix.is_in_strict_zeilenstufenform(), false);
        assert_eq!(matrix.is_in_exact_solution_form(), false);

        let mut matrix_clone = matrix.clone();
        matrix_clone.turn_into_zeilenstufenform();
        println!("{matrix_clone:#?}");
        matrix_clone.turn_into_strict_zeilenstufenform();
        println!("{matrix_clone:#?}");

        let solution = matrix_clone.solve_as_system_of_linear_equations();
        println!("{solution:#?}");
    }

    #[test]
    fn test_7() {
        let matrix = Matrix::from_cells([[1.0, 2.0, 3.0, 4.0]]);

        let mut matrix_clone = matrix.clone();

        let solution = matrix_clone.solve_as_system_of_linear_equations();
        println!("{solution:#?}");
    }

    #[test]
    fn test_8() {
        let matrix = Matrix::from_cells([[1.0, 2.0, 3.0]])
            * (Matrix::from_cells([[4.0], [0.0], [0.0]])
                + Matrix::from_cells([
                    [1.0, 0.0],
                    [0.0, 0.5],
                    [-0.3333333333333333, -0.3333333333333333],
                ]) * Matrix::from_cells([[-123.0], [423.0]]));

        println!("{matrix:?}");

        /*let mut matrix_clone = matrix.clone();

        let solution = matrix_clone.solve_as_system_of_linear_equations();
        println!("{solution:#?}");*/
    }
}
