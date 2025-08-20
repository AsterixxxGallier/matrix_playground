#![allow(dead_code)]

use std::iter::repeat_with;
use std::ops::{Add, AddAssign, Index, IndexMut};

#[derive(Clone, Eq, PartialEq, Debug)]
pub struct Matrix<T> {
    pub rows: usize,
    pub columns: usize,
    pub cells: Vec<Vec<T>>,
}

impl<T> Matrix<T> {
    pub fn new(rows: usize, columns: usize) -> Self
    where
        T: Default,
    {
        Self {
            rows,
            columns,
            cells: (0..rows)
                .map(|_| repeat_with(T::default).take(columns).collect())
                .collect(),
        }
    }

    pub fn from_cells<const ROWS: usize, const COLUMNS: usize>(
        cells: [[T; COLUMNS]; ROWS],
    ) -> Self {
        Self {
            rows: ROWS,
            columns: COLUMNS,
            cells: cells
                .into_iter()
                .map(|row| {
                    let mut row_vec = Vec::with_capacity(COLUMNS);
                    for cell in row {
                        row_vec.push(cell);
                    }
                    row_vec
                })
                .collect(),
        }
    }

    pub fn add_rows(&mut self, rows: usize) where T: Default {
        self.rows += rows;
        self.cells.extend((0..rows).map(|_| (0..self.columns).map(|_| T::default()).collect()));
    }

    pub fn add_columns(&mut self, columns: usize) where T: Default {
        self.columns += columns;
        for row_cells in &mut self.cells {
            row_cells.extend((0..columns).map(|_| T::default()));
        }
    }

    pub fn attach_below(&mut self, mut other: Matrix<T>) {
        assert_eq!(self.columns, other.columns);
        self.rows += other.rows;
        self.cells.append(&mut other.cells);
    }

    pub fn attach_to_the_right(&mut self, other: Matrix<T>) {
        assert_eq!(self.rows, other.rows);
        self.columns += other.columns;
        for (this_row, mut other_row) in self.cells.iter_mut().zip(other.cells) {
            this_row.append(&mut other_row);
        }
    }

    pub fn row(&self, row: usize) -> Matrix<T>
    where
        T: Clone,
    {
        Self {
            rows: 1,
            columns: self.columns,
            cells: vec![self.cells[row].clone()],
        }
    }

    pub fn column(&self, column: usize) -> Matrix<T>
    where
        T: Clone,
    {
        Self {
            rows: self.rows,
            columns: 1,
            cells: self
                .cells
                .iter()
                .map(|row| vec![row[column].clone()])
                .collect(),
        }
    }

    pub fn transpose(&mut self)
    where
        T: Clone,
    {
        let rows = self.columns;
        let columns = self.rows;
        let old_cells = &self.cells;
        let cells = (0..rows)
            .map(move |row| {
                (0..columns)
                    .map(move |column| old_cells[column][row].clone())
                    .collect()
            })
            .collect();

        self.rows = rows;
        self.columns = columns;
        self.cells = cells;
    }
    
    pub fn swap_rows(&mut self, row_one: usize, row_two: usize) {
        self.cells.swap(row_one, row_two);
    }
    
    pub fn swap_columns(&mut self, column_one: usize, column_two: usize) {
        for row_cells in &mut self.cells {
            row_cells.swap(column_one, column_two);
        }
    }
}

impl<T> Index<(usize, usize)> for Matrix<T> {
    type Output = T;

    fn index(&self, (row, column): (usize, usize)) -> &Self::Output {
        &self.cells[row][column]
    }
}

impl<T> IndexMut<(usize, usize)> for Matrix<T> {
    fn index_mut(&mut self, (row, column): (usize, usize)) -> &mut Self::Output {
        &mut self.cells[row][column]
    }
}

impl<T> Add for Matrix<T> where T: AddAssign {
    type Output = Self;

    fn add(mut self, other: Self) -> Self::Output {
        assert_eq!(self.rows, other.rows);
        assert_eq!(self.columns, other.columns);
        for (this_rows, other_rows) in self.cells.iter_mut().zip(other.cells) {
            for (this_cell, other_cell) in this_rows.iter_mut().zip(other_rows) {
                *this_cell += other_cell;
            }
        }
        self
    }
}